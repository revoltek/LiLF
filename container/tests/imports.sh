#!/usr/bin/env bash

# Inside the container, clone LiLF, cd to its directory and paste the following:

python - <<'PY'
import ast
import pathlib
import importlib.util
import sys

root = pathlib.Path(".")

stdlib = set(sys.stdlib_module_names)

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

            if importlib.util.find_spec(topmod) is None:

                print(f"\nMISSING: {topmod}")
                print(f"FILE:    {path}")
                print(f"LINE:    {node.lineno}")

                line = txt.splitlines()[node.lineno - 1]
                print(f"CODE:    {line.strip()}")
PY
