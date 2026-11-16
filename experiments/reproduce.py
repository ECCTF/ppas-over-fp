#!/usr/bin/env python3
"""Rebuild and independently audit the complete p=47,71,127 rational graphs.

Run from any directory with Sage's Python, for example::

    sage -python experiments/reproduce.py --primes 47 71 127

Results replace the selected primes' files in experiments/results. Compilation
products, the execution log, and a completion summary are written to build/.
Every subprocess must succeed before the next stage starts; in particular, the
independent exhaustive audit precedes the Laplacian check and visual exports.
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time


BASE = Path(__file__).resolve().parent
ROOT = BASE.parent
OUT = BASE / "results"
BUILD = ROOT / "build"
SUPPORTED_PRIMES = (47, 71, 127)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--primes", nargs="+", type=int, choices=SUPPORTED_PRIMES,
                        default=list(SUPPORTED_PRIMES),
                        help="published examples to rebuild (default: 47 71 127)")
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"),
                        help="C++17 compiler executable or path (default: CXX, then g++)")
    return parser.parse_args()


def run(command: list[str], log) -> None:
    """Stream a command's combined output to the terminal and persistent log."""
    banner = "\n$ " + shlex.join(command) + "\n"
    print(banner, end="", flush=True)
    log.write(banner)
    log.flush()
    started = time.monotonic()
    with subprocess.Popen(command, cwd=BASE, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                          errors="replace", bufsize=1) as process:
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        returncode = process.wait()
    elapsed = f"[exit {returncode}; {time.monotonic() - started:.1f} seconds]\n"
    print(elapsed, end="", flush=True)
    log.write(elapsed)
    log.flush()
    if returncode:
        raise subprocess.CalledProcessError(returncode, command)


def main() -> int:
    args = parse_args()
    if sys.version_info < (3, 10):
        raise SystemExit("Python 3.10 or later is required.")
    if sys.flags.optimize:
        raise SystemExit("Do not use Python -O or PYTHONOPTIMIZE: verification uses assertions.")
    try:
        importlib.import_module("sage.all")
    except ImportError as error:
        raise SystemExit("Run this script with Sage: sage -python experiments/reproduce.py") from error
    for module in ("numpy", "matplotlib", "networkx"):
        try:
            importlib.import_module(module)
        except ImportError as error:
            raise SystemExit(f"Missing dependency {module!r}; install the repository environment.") from error
    compiler = shutil.which(args.cxx)
    if compiler is None:
        raise SystemExit(f"C++ compiler not found: {args.cxx!r}. Set CXX or pass --cxx PATH.")

    primes = list(dict.fromkeys(args.primes))
    OUT.mkdir(parents=True, exist_ok=True)
    BUILD.mkdir(parents=True, exist_ok=True)
    executable = BUILD / ("ppas_branch_census.exe" if os.name == "nt" else "ppas_branch_census")
    log_path = BUILD / "reproduce.log"
    summary_path = BUILD / "reproduction_summary.json"
    summary = {"primes": primes, "complete": False, "completed_primes": []}
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    started = time.monotonic()

    def python(script: str, *arguments: object) -> list[str]:
        return [sys.executable, "-u", str(BASE / script), *map(str, arguments)]

    try:
        with log_path.open("w", encoding="utf-8") as log:
            run([compiler, "-O3", "-std=c++17", str(BASE / "enumerate_branch_census.cpp"),
                 "-o", str(executable)], log)
            for p in primes:
                for command in (
                    python("bootstrap_explicit_graphs.py", p),
                    python("split_gluing_seeds.py", p),
                    python("restriction_gluing_seeds.py", p),
                    python("rational_census.py", p),
                    python("assemble_rational_census.py", p),
                    python("verify_gluing_seeds.py", p),
                    python("verify_rational_census.py", OUT / f"p{p}_rational_census.json"),
                    [str(executable), str(p), str(OUT / f"p{p}_exhaustive_branch_models.csv")],
                    python("audit_explicit_graphs.py", p),
                    python("verify_explicit_laplacians.py", p),
                    python("export_explicit_graphs.py", p),
                ):
                    run(command, log)
                audit = json.loads((OUT / f"p{p}_full_enumeration_audit.json").read_text(encoding="utf-8"))
                rank = json.loads((OUT / f"p{p}_laplacian_check.json").read_text(encoding="utf-8"))
                summary["completed_primes"].append({
                    "p": p, "vertices": audit["vertices"],
                    "directed_kernels": audit["directed_kernels"],
                    "bfs_components": audit["bfs_components"],
                    "laplacian_rank_Q": rank["rank_Q"],
                })
                summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    except subprocess.CalledProcessError as error:
        print(f"\nReproduction stopped after a failed stage. See {log_path}", file=sys.stderr)
        return error.returncode if error.returncode > 0 else 1

    summary.update(complete=True, elapsed_seconds=round(time.monotonic() - started, 3))
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\nCompleted all requested graphs. Results: {OUT}\nLog: {log_path}\nSummary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
