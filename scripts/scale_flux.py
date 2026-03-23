#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


# -------------------------
# Helpers
# -------------------------
def run(cmd: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    print("\n" + "=" * 80)
    print(" ".join(cmd))
    if cwd is not None:
        print(f"(cwd: {cwd})")
    subprocess.run(cmd, check=True, cwd=None if cwd is None else cwd.as_posix(), env=env)


def list_tc_ms_dirs(ms_root: Path) -> list[Path]:
    if not ms_root.exists():
        raise FileNotFoundError(f"Expected folder '{ms_root}/' in the current working directory.")

    ms_paths = sorted(p for p in ms_root.glob("TC0*.MS") if p.is_dir())
    if not ms_paths:
        raise FileNotFoundError(f"No MeasurementSets found matching {ms_root}/TC0*.MS")

    out: list[Path] = []
    for p in ms_paths:
        if re.match(r"^TC0\d\.MS$", p.name):
            out.append(p)
    return out


def infer_c_from_model(model_path: Path) -> str:
    """
    Infer 'c0' or 'c1' from model filename like tgts-c0.skymodel / tgts-c1.skymodel.
    """
    m = re.search(r"-c([01])(\.|$)", model_path.name)
    if not m:
        raise ValueError(f"Could not infer c0/c1 from model filename: {model_path.name}")
    return f"c{m.group(1)}"


# -------------------------
# Skymodel patch parsing
# -------------------------
def split_fields(line: str) -> list[str]:
    return [p.strip() for p in line.rstrip("\n").split(",")]


def is_patch_definition_line(line: str) -> bool:
    parts = split_fields(line)
    return len(parts) == 5 and parts[0] == "" and parts[1] == "" and parts[2] != ""


def read_patch_names_from_skymodel(skymodel_path: Path) -> list[str]:
    text = skymodel_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    patches: list[str] = []
    for ln in lines:
        if is_patch_definition_line(ln):
            parts = split_fields(ln)
            patches.append(parts[2])

    patches = [p for p in patches if re.fullmatch(r"patch_\d\d", p) is not None]

    seen = set()
    uniq: list[str] = []
    for p in patches:
        if p not in seen:
            seen.add(p)
            uniq.append(p)

    if not uniq:
        raise RuntimeError(
            f"No patch_XX patch definition lines found in {skymodel_path}.\n"
            "If your skymodel does not have explicit patch definition rows, the parser must be adapted."
        )
    return uniq


# -------------------------
# TAQL: sum per-patch columns -> MODEL_DATA
# -------------------------
def taql_sum_model_columns(ms_path: Path, patch_cols: list[str], out_col: str = "MODEL_DATA") -> None:
    if not patch_cols:
        raise ValueError("No patch columns provided for TAQL sum.")
    expr = " + ".join(patch_cols)
    query = f"UPDATE {ms_path.as_posix()} SET {out_col} = {expr}"
    run(["taql", query])


# -------------------------
# Merge + LoSoTo once (global scaling)
# -------------------------
def merge_and_run_losoto(
    *,
    h5s: list[Path],
    out_h5: Path,
    losoto_parsets: list[Path],
    plots_dir: Path,
    solset: str = "sol000",
    workdir: Path = Path("."),
) -> Path:
    if not h5s:
        raise RuntimeError("No H5 files provided to merge.")

    workdir = workdir.resolve()
    out_h5 = out_h5.resolve()
    plots_dir = plots_dir.resolve()

    # Merge (or copy) into out_h5
    if len(h5s) > 1:
        # remove existing output (file or directory)
        if out_h5.exists():
            if out_h5.is_dir():
                shutil.rmtree(out_h5)
            else:
                out_h5.unlink()

        run(
            [
                "H5parm_collector.py",
                "-V",
                "-s",
                solset,
                "-o",
                out_h5.as_posix(),
                *[p.resolve().as_posix() for p in h5s],
            ],
            cwd=workdir,
        )
    else:
        src = h5s[0].resolve()
        if out_h5.exists():
            if out_h5.is_dir():
                shutil.rmtree(out_h5)
            else:
                out_h5.unlink()
        # h5parm is usually a directory; copytree handles that
        shutil.copytree(src, out_h5)

    # Clean plots dir and run LoSoTo (writes into workdir/plots)
    tmp_plots = workdir / "plots"
    if tmp_plots.exists():
        shutil.rmtree(tmp_plots)
    tmp_plots.mkdir(parents=True, exist_ok=True)

    for parset in losoto_parsets:
        run(["losoto", "-V", out_h5.as_posix(), parset.resolve().as_posix()], cwd=workdir)

    # Move plots/* to plots_dir
    if plots_dir.exists():
        shutil.rmtree(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    if tmp_plots.exists():
        for item in tmp_plots.iterdir():
            shutil.move(item.as_posix(), (plots_dir / item.name).as_posix())
        shutil.rmtree(tmp_plots)

    return out_h5


# -------------------------
# Per-MS pipeline (runs in its own work dir)
# Now returns the per-MS solved h5parm path, but does NOT run LoSoTo or apply amplitudeSmooth.
# -------------------------
def pipeline_one_ms(
    *,
    ms_path: Path,
    skymodel: Path,
    patch_names: list[str],
    parset_predict: Path,
    parset_cor: Path,
    parset_solve: Path,
    cal_fr_h5: Path,
    cal_tec_h5: Path,
    antenna_constraint_value: str,
    solint: int,
    work_root: Path,
    dp3_threads_env: dict[str, str],
    do_fr_corruption: bool,
    solve_datacolumn: str,  # "DATA" or "CORRECTED_DATA"
) -> Path:
    m = re.match(r"^(TC0\d)\.MS$", ms_path.name)
    if not m:
        raise ValueError(f"Unexpected MS name: {ms_path.name}")
    tcid = m.group(1)

    # Make per-MS working directory to avoid collisions
    workdir = work_root / f"work-{tcid}"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    # Absolute paths so cwd doesn't matter
    ms_path = ms_path.resolve()
    skymodel = skymodel.resolve()
    parset_predict = parset_predict.resolve()
    parset_cor = parset_cor.resolve()
    parset_solve = parset_solve.resolve()
    cal_fr_h5 = cal_fr_h5.resolve()
    cal_tec_h5 = cal_tec_h5.resolve()

    patch_cols = patch_names  # patch_00 etc

    # --- Per direction: predict -> (optional FR/RM) -> TEC phase
    for patch in patch_names:
        # Predict beam into column "patch_XX"
        run(
            [
                "DP3",
                parset_predict.as_posix(),
                f"msin={ms_path.as_posix()}",
                f"pre.sourcedb={skymodel.as_posix()}",
                f"pre.sources={patch}",
                f"msout.datacolumn={patch}",
                "pre.correctfreqsmearing=True",
            ],
            env=dp3_threads_env,
        )

        # FR / RM corruption only for c0
        if do_fr_corruption:
            run(
                [
                    "DP3",
                    parset_cor.as_posix(),
                    f"msin={ms_path.as_posix()}",
                    f"cor.parmdb={cal_fr_h5.as_posix()}",
                    f"msin.datacolumn={patch}",
                    f"msout.datacolumn={patch}",
                    "cor.correction=rotationmeasure000",
                    "cor.invert=False",
                ],
                env=dp3_threads_env,
            )

        # TEC phase corruption
        run(
            [
                "DP3",
                parset_cor.as_posix(),
                f"msin={ms_path.as_posix()}",
                f"msin.datacolumn={patch}",
                f"msout.datacolumn={patch}",
                f"cor.direction=[{patch}]",
                f"cor.parmdb={cal_tec_h5.as_posix()}",
                "cor.correction=phase000",
                "cor.invert=False",
            ],
            env=dp3_threads_env,
        )

    # --- Sum -> MODEL_DATA
    taql_sum_model_columns(ms_path, patch_cols, out_col="MODEL_DATA")

    # --- Solve amplitudes (per-MS)
    h5parm = (workdir / f"amp-di_{tcid}.h5").resolve()
    run(
        [
            "DP3",
            parset_solve.as_posix(),
            f"msin={ms_path.as_posix()}",
            f"msin.datacolumn={solve_datacolumn}",
            "sol.datause=full",
            "sol.nchan=16",
            "sol.modeldatacolumns=[MODEL_DATA]",
            "sol.mode=diagonal",
            f"sol.h5parm={h5parm.as_posix()}",
            f"sol.solint={solint}",
            f"sol.antennaconstraint={antenna_constraint_value}",
        ],
        env=dp3_threads_env,
    )

    return h5parm


# -------------------------
# Imaging (serial, after all MS done)
# -------------------------
def run_wsclean(ms_paths: list[Path]) -> None:
    Path("img_amp").mkdir(parents=True, exist_ok=True)

    cmd = [
        "wsclean",
        "-j",
        "256",
        "-beam-size",
        "45asec",
        "-reorder",
        "-parallel-reordering",
        "4",
        "-name",
        "img_amp/wideDDP-c1_amp",
        "-no-update-model-required",
        "-nmiter",
        "12",
        "-auto-threshold",
        "2.0",
        "-auto-mask",
        "4.0",
        "-apply-facet-beam",
        "-facet-beam-update",
        "120",
        "-use-differential-lofar-beam",
        "-local-rms",
        "-local-rms-window",
        "50",
        "-local-rms-strength",
        "0.5",
        "-data-column",
        "CORRECTED_DATA_AMPL",
        "-size",
        "10000",
        "10000",
        "-scale",
        "9arcsec",
        "-weight",
        "briggs",
        "-0.5",
        "-niter",
        "1000000",
        "-gridder",
        "wgridder",
        "-parallel-gridding",
        "4",
        "-minuv-l",
        "30",
        "-mgain",
        "0.85",
        "-parallel-deconvolution",
        "32",
        "-join-channels",
        "-fit-spectral-pol",
        "3",
        "-channels-out",
        "4",
        "-deconvolution-channels",
        "3",
        "-multiscale",
        "-multiscale-scale-bias",
        "0.65",
        "-pol",
        "i",
        "-facet-regions",
        "ddparallel/solutions/facetsP-c1.reg",
        "-apply-facet-solutions",
        "ddparallel/solutions/cal-tec-c1.h5",
        "phase000",
    ]
    cmd.extend([p.resolve().as_posix() for p in ms_paths])
    run(cmd)


# -------------------------
# Main
# -------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("model", type=Path, help="Input skymodel (e.g. tgts-c0.skymodel or tgts-c1.skymodel)")
    ap.add_argument("--ms-root", type=Path, default=Path("mss"))
    ap.add_argument("--workers", type=int, default=2, help="Number of MS to process in parallel")
    ap.add_argument("--solint", type=int, default=75)
    ap.add_argument("--work-root", type=Path, default=Path("."), help="Where to create work-TC0n dirs")
    args = ap.parse_args()

    c_tag = infer_c_from_model(args.model)  # 'c0' or 'c1'

    # Rules you requested:
    # - TEC solution follows c_tag (cal-tec-c0 or cal-tec-c1)
    # - If c1: do NOT corrupt model with FR, and solve amplitudes using CORRECTED_DATA (not DATA)
    do_fr_corruption = c_tag == "c0"
    solve_datacolumn = "DATA" if c_tag == "c0" else "CORRECTED_DATA"

    ms_paths = list_tc_ms_dirs(args.ms_root)

    patch_names = read_patch_names_from_skymodel(args.model)
    print(f"Model tag: {c_tag}")
    print(f"FR corruption: {do_fr_corruption}")
    print(f"Amplitude solve datacolumn: {solve_datacolumn}")
    print(f"Found {len(patch_names)} patches in model: {patch_names}")

    parset_predict = Path("pipeline-old/parsets/LOFAR_ddparallel/DP3-predict-beam.parset")
    parset_cor = Path("pipeline-old/parsets/LOFAR_ddparallel/DP3-cor.parset")
    parset_solve = Path("pipeline-old/parsets/LOFAR_ddparallel/DP3-soldd.parset")

    losoto_plot_amp = Path("pipeline-old/parsets/LOFAR_ddparallel/losoto-plot-amp.parset")
    losoto_plot_ph = Path("pipeline-old/parsets/LOFAR_ddparallel/losoto-plot-ph.parset")
    losoto_amp_di = Path("pipeline-old/parsets/LOFAR_ddparallel/losoto-amp-di.parset")

    cal_fr_h5 = Path("ddparallel/solutions/cal-fr.h5")
    cal_tec_h5 = Path(f"ddparallel/solutions/cal-tec-{c_tag}.h5")

    antenna_constraint_value = (
        "[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,"
        "CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,"
        "CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,"
        "RS106LBA,RS205LBA,RS305LBA,RS306LBA,RS503LB]]"
    )

    # Avoid oversubscription when running multiple DP3 processes
    base_env = dict(os.environ)
    base_env["OMP_NUM_THREADS"] =  "1"
    base_env["OPENBLAS_NUM_THREADS"] = "1"
    base_env["MKL_NUM_THREADS"] = "1"
    base_env["NUMEXPR_NUM_THREADS"] = "1"

    # Parallel per-MS: solve only, collect per-MS h5parms
    solved_h5s: list[Path] = []
    failures: list[tuple[Path, str]] = []

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        fut_to_ms: dict = {}
        for ms in ms_paths:
            fut = ex.submit(
                pipeline_one_ms,
                ms_path=ms,
                skymodel=args.model,
                patch_names=patch_names,
                parset_predict=parset_predict,
                parset_cor=parset_cor,
                parset_solve=parset_solve,
                cal_fr_h5=cal_fr_h5,
                cal_tec_h5=cal_tec_h5,
                antenna_constraint_value=antenna_constraint_value,
                solint=args.solint,
                work_root=args.work_root,
                dp3_threads_env=base_env,
                do_fr_corruption=do_fr_corruption,
                solve_datacolumn=solve_datacolumn,
            )
            fut_to_ms[fut] = ms

        for fut in as_completed(fut_to_ms):
            ms = fut_to_ms[fut]
            try:
                h5 = fut.result()
                solved_h5s.append(Path(h5))
            except Exception as e:
                failures.append((ms, repr(e)))

    if failures:
        print("\nSome MS jobs failed:")
        for ms, err in failures:
            print(f"  {ms}: {err}")
        raise SystemExit(1)

    # Merge + run LoSoTo 
    merged_h5 = merge_and_run_losoto(
        h5s=sorted(solved_h5s),
        out_h5=Path(f"cal-amp-di-{c_tag}.h5"),
        losoto_parsets=[losoto_plot_amp, losoto_plot_ph, losoto_amp_di],
        plots_dir=Path(f"plots-amp-di-{c_tag}"),
        solset="sol000",
        workdir=Path("."),
    )
    print(f"\nMerged+scaled h5parm: {merged_h5}")

    # Apply amplitudeSmooth using the merged/scaled h5 to each MS
    for ms in ms_paths:
        run(
            [
                "DP3",
                parset_cor.as_posix(),
                f"msin={ms.resolve().as_posix()}",
                "msin.datacolumn=CORRECTED_DATA",
                "msout.datacolumn=CORRECTED_DATA_AMPL",
                f"cor.parmdb={merged_h5.as_posix()}",
                "cor.correction=amplitudeSmooth",
                "cor.updateweights=False",
            ],
            env=base_env,
        )

    # Imaging after all MS done
    run_wsclean(ms_paths)


if __name__ == "__main__":
    main()