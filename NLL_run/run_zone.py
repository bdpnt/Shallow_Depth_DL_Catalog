"""
run_zone.py
============================
Run the full NonLinLoc pipeline (Vel2Grid → Grid2Time → NLLoc) for a single
zone given a .in run file.

For a first-pass run file (run_<N>_DELAYS.in) all three programs are executed
in sequence.  For a corrections-pass file (run_<N>_NLL.in) only NLLoc is run
because the velocity and travel-time grids were already computed in the first
pass.

Usage
-----
    python NLL_run/run_zone.py run/nll/run_1_DELAYS.in
    python NLL_run/run_zone.py run/nll/run_1_NLL.in
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
# NLL binary directory
# ---------------------------------------------------------------------------

# TEMPORARY: locally patched build ('*' picks in EDT_OT_WT), until the fix is
# upstream -- see NonLinLoc_fixed/README_TEMPORARY.md. Revert to NonLinLoc/src/bin.
_NLL_BIN = "/Users/bdupont/Desktop/Codes/NonLinLoc_fixed/src/bin"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

logging.basicConfig(
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
    level=logging.INFO,
)
log = logging.getLogger(__name__)

_MIN_FREE_GB = 5.0

_FATAL_PATTERNS = (
    "no space left on device",
    "cannot allocate memory",
    "out of memory",
    "malloc failed",
    "malloc: failed",
)


def _exe(name: str) -> str:
    return os.path.join(_NLL_BIN, name)


def _check_disk(run_in: str) -> None:
    out_dir = os.path.dirname(run_in)
    with open(run_in) as f:
        for line in f:
            if line.startswith("VGOUT"):
                parts = line.split()
                if len(parts) >= 2:
                    out_dir = os.path.dirname(parts[1]) or "."
                break
    try:
        free_gb = shutil.disk_usage(out_dir).free / 1e9
        if free_gb < _MIN_FREE_GB:
            log.warning("Low disk space: %.1f GB free at %s", free_gb, out_dir)
    except OSError:
        pass


def _run(cmd: list[str], label: str, zone_label: str = "") -> None:
    prefix = f"[{zone_label}] " if zone_label else ""
    log.info("%sStarting %s", prefix, label)
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )
    lines: list[str] = []
    for line in proc.stdout:
        sys.stdout.write(f"{prefix}{line}")
        sys.stdout.flush()
        lines.append(line)
    proc.wait()

    if proc.returncode != 0:
        log.error("%s%s failed (exit %d)", prefix, label, proc.returncode)
        sys.exit(proc.returncode)

    output_lower = "".join(lines).lower()
    for pattern in _FATAL_PATTERNS:
        if pattern in output_lower:
            log.error("%s%s: fatal error detected in output (%r)", prefix, label, pattern)
            sys.exit(1)

    log.info("%s%s done", prefix, label)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_zone(run_in: str, *, corrections_pass: bool = False, zone_label: str = "") -> None:
    """Run Vel2Grid → Grid2Time → NLLoc for the zone described by *run_in*.

    Set *corrections_pass* to True to skip Vel2Grid and Grid2Time (grids
    already built during the first pass). *zone_label* prefixes log output,
    useful when multiple zones run concurrently.
    """
    run_in = os.path.abspath(run_in)
    if not os.path.isfile(run_in):
        raise FileNotFoundError(f"Run file not found: {run_in}")

    if not corrections_pass:
        _check_disk(run_in)
        _run([_exe("Vel2Grid"),   run_in], "Vel2Grid",  zone_label)
        _run([_exe("Grid2Time"),  run_in], "Grid2Time", zone_label)

    _run([_exe("NLLoc"), run_in], "NLLoc", zone_label)


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
