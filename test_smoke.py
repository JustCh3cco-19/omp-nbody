#!/usr/bin/env python3
"""Local correctness smoke tests."""
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
ENV = os.environ.copy()
ENV.setdefault("OMP_NUM_THREADS", "4")
FULL = os.environ.get("FULL", "0") == "1"


def run(binary, args, cwd, mpi=False):
    launcher = shlex.split(os.environ.get("MPIEXEC", "mpirun -np 4")) if mpi else []
    result = subprocess.run(launcher + [str(ROOT / binary)] + list(map(str, args)),
                            cwd=cwd, env=ENV, text=True, capture_output=True,
                            timeout=float(os.environ.get("TEST_TIMEOUT", "300")))
    if result.returncode:
        raise RuntimeError(f"{binary} exited {result.returncode}\n{result.stderr}")
    return result.stdout


def compare(actual, expected):
    a, b = actual.split(), expected.split()
    if not a or len(a) != len(b):
        raise ValueError(f"output length mismatch: {len(a)} vs {len(b)} tokens")
    for index, (x, y) in enumerate(zip(a, b)):
        x, y = float(x), float(y)
        if not (math.isfinite(x) and math.isfinite(y)
                and math.isclose(x, y, rel_tol=5e-5, abs_tol=1e-9)):
            raise ValueError(f"numeric mismatch at token {index}: {x} vs {y}")


def main():
    with tempfile.TemporaryDirectory(prefix="mcp-smoke-") as temporary:
        work = Path(temporary)
        for bodies in ((20000, 55000) if FULL else (20000,)):
            name = f"positions.xyz.{bodies}"
            (work / name).symlink_to(ROOT / name)
            (work / "results").unlink(missing_ok=True)
            run("nbody_omp.out", [bodies], work)
            compare((work / "results").read_text(), (ROOT / f"results.{bodies}").read_text())
            print(f"PASS: {bodies} bodies")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
