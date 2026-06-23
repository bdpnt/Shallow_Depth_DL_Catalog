"""
run_zone.py
============================
Run the full NonLinLoc pipeline (Vel2Grid → Grid2Time → NLLoc) for a single
zone given a .in run file.

For a first-pass run file (run_<N>.in) all three programs are executed in
sequence.  For a corrections-pass file (run_<N>_PR.in) only NLLoc is run
because the velocity and travel-time grids were already computed in the first
pass.

Usage
-----
    python NLL_run/run_zone.py run/run_1.in
    python NLL_run/run_zone.py run/run_1_PR.in
"""

import argparse
import logging
import os
import subprocess
import sys

# ---------------------------------------------------------------------------
# NLL binary directory
# ---------------------------------------------------------------------------

_NLL_BIN = "/Users/bdupont/Desktop/Codes/NonLinLoc/src/bin"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

logging.basicConfig(
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
    level=logging.INFO,
)
log = logging.getLogger(__name__)


def _exe(name: str) -> str:
    return os.path.join(_NLL_BIN, name)


def _run(cmd: list[str], label: str) -> None:
    log.info("Starting %s", label)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        log.error("%s failed (exit %d)", label, result.returncode)
        sys.exit(result.returncode)
    log.info("%s done", label)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_zone(run_in: str, *, corrections_pass: bool = False) -> None:
    """Run Vel2Grid → Grid2Time → NLLoc for the zone described by *run_in*.

    Set *corrections_pass* to True to skip Vel2Grid and Grid2Time (grids
    already built during the first pass).
    """
    run_in = os.path.abspath(run_in)
    if not os.path.isfile(run_in):
        raise FileNotFoundError(f"Run file not found: {run_in}")

    if not corrections_pass:
        _run([_exe("Vel2Grid"),   run_in], "Vel2Grid")
        _run([_exe("Grid2Time"),  run_in], "Grid2Time")

    _run([_exe("NLLoc"), run_in], "NLLoc")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _main() -> None:
    parser = argparse.ArgumentParser(
        description="Run Vel2Grid → Grid2Time → NLLoc for one zone."
    )
    parser.add_argument("run_in", help="Path to the NLL .in run file")
    parser.add_argument(
        "--corrections-pass",
        action="store_true",
        help="Skip Vel2Grid / Grid2Time (grids already built in first pass)",
    )
    args = parser.parse_args()
    run_zone(args.run_in, corrections_pass=args.corrections_pass)


if __name__ == "__main__":
    _main()
