#!/usr/bin/env bash

# Inside the container, clone LiLF, cd to its directory and paste the following:

python - <<'PY'
import ast
import pathlib
import importlib
import warnings
import sys

# Hide SyntaxWarning: invalid escape sequence '\['
warnings.filterwarnings(
    "ignore",
    message=r"invalid escape sequence",
    category=SyntaxWarning,
)

root = pathlib.Path(".")
stdlib = set(sys.stdlib_module_names)
checked = set()

for path in root.rglob("*.py"):
    try:
        txt = path.read_text(errors="ignore")
        tree = ast.parse(txt)
    except Exception:
        continue

    for node in ast.walk(tree):
        modules = []
        if isinstance(node, ast.Import):
            for n in node.names:
                modules.append((n.name.split(".")[0], n.name))

        elif isinstance(node, ast.ImportFrom):
            if node.module:
                modules.append(
                    (node.module.split(".")[0], node.module)
                )

        for topmod, fullname in modules:
            if topmod in stdlib:
                continue

            key = (topmod, path, node.lineno)
            if key in checked:
                continue
            checked.add(key)
            line = txt.splitlines()[node.lineno - 1].strip()

            try:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always")
                    importlib.import_module(topmod)
                    if w:
                        print(f"\nWARNING: {topmod}")
                        print(f"FILE:    ./{path}")
                        print(f"LINE:    {node.lineno}")
                        print(f"CODE:    {line}")

                        for warn in w:
                            print(
                                f"WARN:    "
                                f"{warn.category.__name__}: "
                                f"{warn.message}"
                            )

            except ModuleNotFoundError:
                print(f"\nMISSING: {topmod}")
                print(f"FILE:    ./{path}")
                print(f"LINE:    {node.lineno}")
                print(f"CODE:    {line}")

            except Exception as e:
                print(f"\nERROR:   {topmod}")
                print(f"FILE:    ./{path}")
                print(f"LINE:    {node.lineno}")
                print(f"CODE:    {line}")
                print(f"EXC:     {type(e).__name__}: {e}")

PY
