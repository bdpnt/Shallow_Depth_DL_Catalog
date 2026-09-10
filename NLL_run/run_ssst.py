"""
run_ssst.py
============================
Iterative NLLoc + Loc2ssst SSST relocation workflow for one zone.

Port of CODES_SSST/run_ssst_VF.py (itself a Python port of run_ssst_VF.bash
by Anthony Lomax, ALomax Scientific). All the former module-level constants
are now fields of the SSSTRunParams dataclass (Python calls) with matching
CLI options (command-line calls) — configuring a run never requires editing
this module.

Iteration N (0-based) = one parallel NLLoc relocation + one parallel Loc2ssst
computing correction grids smoothed with charDists[N] (first value huge ->
static corrections). Iteration len(charDists) is a final NLLoc-only
relocation with SAVE_NLLOC_OCTREE and SAVE_FMAMP forced on. The iteration
window [iteration_start, iteration_stop] allows resuming an interrupted
campaign or running only the first iterations now and the rest later.

Differences from the validated CODES_SSST original (porting requirements of
SSST_INTEGRATION.md §8): fail-fast on any non-zero subprocess exit and on
empty obs/time-grid globs, round-robin obs chunking, separate NLLoc/Loc2ssst
core counts, per-iteration chunk-directory cleanup, per-zone tmp dir /
journal / ssst.list, and binaries invoked from the NLL bin directory instead
of relying on PATH. build_grids() prepares the initial P and S grids
(Vel2Grid + Grid2Time twice, once per phase).

Usage
-----
    python NLL_run/run_ssst.py \\
        --run-name ssst_run1 --project-name Pyrenees_1 --model Pyrenees_1 \\
        --catalog GLOBAL_1 \\
        --control-nll  run/ssst/run_1_NLL.in \\
        --control-ssst run/ssst/run_1_SSST.in \\
        --loc-obs  'obs/nlloc_obs/GLOBAL_1/*.nlloc_obs' \\
        --stations stations/GTSRCE_SSST_1.txt \\
        --init-time-root run/ssst_time/Pyrenees_1/Pyrenees_1 \\
        --out run/ssst_loc --tmp-dir run/ssst/tmp_1 \\
        --build-grids
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

# NLL binary directory (same as run_zone.py)
_NLL_BIN = "/Users/bdupont/Desktop/Codes/NonLinLoc/src/bin"

_MIN_FREE_GB = 50.0

_FATAL_PATTERNS = (
    "no space left on device",
    "cannot allocate memory",
    "out of memory",
    "malloc failed",
    "malloc: failed",
)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SSSTRunParams:
    runName:      str    # run identifier (output subdir under outRoot)
    projectName:  str    # basename of the NLLoc output files (e.g. Pyrenees_1)
    model:        str    # basename of the velocity/time grids (e.g. Pyrenees_1)
    catalog:      str    # obs set used to derive corrections (e.g. GLOBAL_1)
    fileRunNLL:   str    # NLLoc control file (run/ssst/run_<N>_NLL.in)
    fileRunSSST:  str    # Loc2ssst control file (run/ssst/run_<N>_SSST.in)
    locObs:       str    # glob of per-event obs files (obs/nlloc_obs/GLOBAL_<N>/*.nlloc_obs)
    stationsFile: str    # GTSRCE file: station list split across Loc2ssst instances
    initTimeRoot: str    # initial time-grid file root (run/ssst_time/Pyrenees_<N>/Pyrenees_<N>)
    outRoot:      str    # output root (run/ssst_loc)
    tmpDir:       str    # per-zone scratch dir (run/ssst/tmp_<N>)
    logFile:      str    # per-zone run journal (run/ssst/run_ssst_<N>.log)
    ssstListFile: str    # list of final summary-hyp paths (run/ssst/ssst.list)

    catalogFinal:  str = None    # obs set for the final relocation (default: catalog)
    locObsFinal:   str = None    # glob for the final relocation (default: locObs)
    ssstModelName: str = None    # output subdir name (default: <model>_SSST)

    # SSST smoothing schedule (km); len(charDists) SSST iterations + 1 final
    charDists: list  = field(default_factory=lambda: [9999, 50, 15, 5, 1])
    vpvs:      float = -9.99     # LOCMETH VpVsRatio, iteration 0 (-9.99 = real S grids)
    vpvsIter:  float = -9.99     # LOCMETH VpVsRatio once corrections are active

    nllocCores:    int = 10      # parallel NLLoc chunks
    loc2ssstCores: int = 10      # parallel Loc2ssst instances (~2 LSGRID buffers each)

    saveNllocOctree: str = ''    # forced to SAVE_NLLOC_OCTREE on the final iteration
    saveFmamp:       str = ''    # forced to SAVE_FMAMP on the final iteration
    zoneLabel:       str = ''    # prefix for console output when zones are scripted

    def __post_init__(self):
        if self.catalogFinal is None:
            self.catalogFinal = self.catalog
        if self.locObsFinal is None:
            self.locObsFinal = self.locObs
        if self.ssstModelName is None:
            self.ssstModelName = f'{self.model}_SSST'


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _exe(name):
    return os.path.join(_NLL_BIN, name)


def _log_append(params, text):
    """Append text to the per-zone run journal."""
    with open(params.logFile, 'a') as f:
        f.write(text)


def _print(params, text):
    prefix = f'[{params.zoneLabel}] ' if params.zoneLabel else ''
    print(f'{prefix}{text}', flush=True)


def _remove_path(path):
    """Delete a file or a directory tree (dispatches to the right call)."""
    if os.path.isdir(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


def _check_disk(out_dir):
    try:
        free_gb = shutil.disk_usage(out_dir).free / 1e9
        if free_gb < _MIN_FREE_GB:
            print(f'WARNING: low disk space: {free_gb:.1f} GB free at {out_dir}')
    except OSError:
        pass


def _scan_log_for_fatal(log_path):
    """Return the fatal pattern found in a subprocess log, or None."""
    try:
        with open(log_path, 'r', errors='replace') as f:
            content = f.read().lower()
    except OSError:
        return None
    for pattern in _FATAL_PATTERNS:
        if pattern in content:
            return pattern
    return None


def _wait_all(params, procs, label):
    """Wait for every (proc, log_path) pair, then fail-fast on any failure.

    All processes are waited on before raising so no orphan keeps writing
    into the output directories. A non-zero exit status or a fatal disk/
    memory pattern in a subprocess log (NLL programs can report those with
    exit code 0) aborts the iteration.
    """
    failures = []
    for proc, log_path in procs:
        status = proc.wait()
        rel_log = os.path.relpath(log_path, _PROJECT_ROOT) if log_path else '-'
        _print(params, f'Finished: {label} PID={proc.pid} status={status} log={rel_log}')
        if status != 0:
            failures.append(f'{label} PID={proc.pid} exit={status} (log: {log_path})')
        elif log_path:
            pattern = _scan_log_for_fatal(log_path)
            if pattern:
                failures.append(f'{label} PID={proc.pid} fatal {pattern!r} (log: {log_path})')
    if failures:
        raise RuntimeError(f'{label} failed:\n  ' + '\n  '.join(failures))


def _run_checked(params, cmd, label, log_path):
    """Run one foreground command with output to log_path; fail-fast."""
    _print(params, f'Starting {label}')
    with open(log_path, 'w') as log_file:
        # cwd=tmpDir (not _PROJECT_ROOT): all paths in the control files are
        # absolute, and NLL programs drop fixed-name byproducts (e.g.
        # Grid2GMT_SSST.cpt) into the cwd — keep those contained in scratch
        proc = subprocess.Popen(cmd, cwd=params.tmpDir,
                                stdout=log_file, stderr=subprocess.STDOUT)
    _wait_all(params, [(proc, log_path)], label)


def _build_skeleton_conf(params):
    """Build a reusable "skeleton" NLLoc control file from fileRunNLL.

    The INCLUDE, LOCFILES and LOCMETH statements are neutralised (commented
    out) because they change on every iteration and every parallel chunk;
    _run_nlloc_parallel() appends fresh, iteration-specific versions to a
    copy of this skeleton for each NLLoc instance. Everything else (TRANS,
    LOCGRID, LOCSEARCH, LOCQUAL2ERR, ...) passes through unchanged.
    """
    with open(params.fileRunNLL) as f:
        lines = f.readlines()
    out = []
    for line in lines:
        # collapse any existing leading '#'s into a single one, so
        # already-commented lines stay commented
        if 'INCLUDE' in line or 'LOCFILES' in line or 'LOCMETH' in line:
            out.append(re.sub(r'^#*', '#', line))
        else:
            out.append(line)
    skeleton_conf = os.path.join(params.tmpDir, 'NLL_SKELETON_SSST.conf')
    with open(skeleton_conf, 'w') as f:
        f.writelines(out)
    return skeleton_conf


def _read_min_num_phases(path):
    """Read the min-phases threshold from the Loc2ssst control file.

    Returns the NRdgsMin value (2nd number) of the LSPHSTAT statement. It is
    reused as the LOCMETH minNumberPhases of the NLLoc steps, so NLLoc never
    locates events that Loc2ssst would reject on phase count (the two
    thresholds cannot drift apart).
    """
    with open(path) as f:
        for line in f:
            fields = line.split()
            if fields and fields[0] == 'LSPHSTAT':
                return int(float(fields[2]))
    raise ValueError(f'no LSPHSTAT statement found in {path}')


def _read_station_codes(path):
    """Read the Loc2ssst station list from a GTSRCE-style station file.

    One station per line, code in the second column. Blank lines and '#'
    comment lines are skipped; duplicate codes are dropped (keeping the
    first occurrence) so the round-robin split stays disjoint.
    """
    codes = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) < 2 or fields[0].startswith('#'):
                continue
            if fields[1] not in codes:
                codes.append(fields[1])
    if not codes:
        raise ValueError(f'no station codes found in {path}')
    return codes


def _run_nlloc_parallel(params, control_file_tmp, out_root, loc_obs, time_root,
                        min_num_phases_loc, vpvs, save_nlloc_octree,
                        save_fmamp, iteration):
    """Locate all events of one iteration with nllocCores parallel NLLoc
    processes, merge their output into out_root, then delete the chunk dirs.
    """
    # --- remove leftover chunk obs directories from a previous iteration/run
    # so a chunk never silently relocates last iteration's obs files
    for stale in glob.glob(os.path.join(params.tmpDir, 'obsfiles_*')):
        _remove_path(stale)

    # --- fail-fast on an empty obs glob: NLLoc chunks would "succeed" on
    # zero events and Loc2ssst would fail later with no matching .hyp files
    obs_files = sorted(glob.glob(loc_obs))
    if not obs_files:
        raise FileNotFoundError(f'no obs files match: {loc_obs}')

    # --- split the obs files round-robin into nllocCores non-empty chunks
    n = params.nllocCores
    chunks = [obs_files[i::n] for i in range(n)]
    chunks = [c for c in chunks if c]

    # --- launch one NLLoc process per chunk (all run concurrently)
    procs = []
    for index, chunk in enumerate(chunks):
        _print(params, f'Running: NLLoc chunk {index} ({len(chunk)} obs files)')

        # per-chunk control file = shared iteration skeleton + statements below
        control_file_idx = f'{control_file_tmp}_{index}'
        shutil.copy2(control_file_tmp, control_file_idx)

        # copy this chunk's obs files into tmpDir/obsfiles_<index>/ so the
        # LOCFILES wildcard of this chunk sees only its own events
        obsfiles_dir = os.path.join(params.tmpDir, f'obsfiles_{index}')
        os.makedirs(obsfiles_dir, exist_ok=True)
        for ofile in chunk:
            shutil.copy2(ofile, obsfiles_dir)

        # append the iteration/chunk-specific statements that were commented
        # out of the skeleton:
        #   LOCCOM    - free-text comment recorded in the output headers
        #   LOCFILES  - obs input pattern, input travel-time grid root
        #               (initial grids at iteration 0, previous iteration's
        #               SSST-corrected grids afterwards), chunk output root
        #   LOCMETH   - EDT location method: min phases and Vp/Vs ratio
        #               (-9.99 = use S travel-time grids, no fixed ratio)
        #   LOCHYPOUT - oct-tree grids and fmamp only on the final iteration;
        #               SAVE_NLLOC_EXPECTATION on every iteration, so the reported
        #               hypocenter is the PDF expectation the confidence ellipsoid
        #               is actually built around, and the residuals Loc2ssst turns
        #               into corrections are the ones evaluated there
        out_root_idx = f'{out_root}_{index}'
        with open(control_file_idx, 'a') as f:
            f.write(f'LOCCOM {params.runName} {params.ssstModelName} '
                    f'loc_ssst_corr{iteration}\n')
            f.write(f'LOCFILES {obsfiles_dir}/*.nlloc_obs NLLOC_OBS  {time_root}  '
                    f'{out_root_idx}/{params.projectName}  0\n')
            f.write(f'LOCMETH EDT_OT_WT 9999.0 {min_num_phases_loc} -1 -1 {vpvs} 145 -1 1\n')
            f.write(f'LOCHYPOUT SAVE_NLLOC_ALL  NLL_FORMAT_VER_2  SAVE_NLLOC_EXPECTATION  '
                    f'{save_nlloc_octree} {save_fmamp}\n')

        # start NLLoc detached; its console output goes to a per-chunk log
        # so the parallel processes do not interleave on the console
        os.makedirs(out_root_idx, exist_ok=True)
        log_path = f'{control_file_idx}.log'
        log_file = open(log_path, 'w')
        proc = subprocess.Popen([_exe('NLLoc'), control_file_idx],
                                cwd=params.tmpDir,
                                stdout=log_file, stderr=subprocess.STDOUT)
        procs.append((proc, log_path, log_file))

    # --- wait for every chunk; fail-fast before any merge
    try:
        _wait_all(params, [(p, lp) for p, lp, _ in procs], 'NLLoc')
    finally:
        for _, _, log_file in procs:
            log_file.close()

    # --- merge step 1: copy every chunk directory into out_root. Per-event
    # files have unique names so they coexist; the per-chunk SUMMARY files
    # share the same name and overwrite each other here - rebuilt below
    os.makedirs(out_root, exist_ok=True)
    chunk_dirs = sorted(glob.glob(f'{out_root}_[0-9]*'))
    for out_root_idx in chunk_dirs:
        shutil.copytree(out_root_idx, out_root, dirs_exist_ok=True)

    # --- merge step 2: rebuild the summary files by explicit concatenation.
    # .hyp and .stations are one self-contained block per event (plain
    # concat); the CSV keeps the header line from the first chunk only
    # (without this, the copytree left only the last chunk's CSV)
    proj = params.projectName
    for suffix in ('sum.grid0.loc.hyp', 'sum.grid0.loc.stations'):
        with open(os.path.join(out_root, f'{proj}.{suffix}'), 'w') as dst:
            for path in sorted(glob.glob(f'{out_root}_[0-9]*/{proj}.{suffix}')):
                with open(path) as src:
                    dst.write(src.read())

    with open(os.path.join(out_root, f'{proj}.sum.grid0.loc.csv'), 'w') as dst:
        for k, path in enumerate(sorted(glob.glob(f'{out_root}_[0-9]*/{proj}.sum.grid0.loc.csv'))):
            with open(path) as src:
                lines = src.readlines()
            dst.writelines(lines if k == 0 else lines[1:])

    # --- leave a marker explaining that the summaries are concatenations
    warning = os.path.join(out_root,
                           'WARNING.concatenated_output_of_multiple_NLLoc_runs.WARNING')
    with open(warning, 'w') as f:
        f.write('\n')

    # --- per-iteration cleanup: the chunk dirs are pure duplicates once
    # merged ([0-9] never matches the merged {catalog} dir itself)
    for out_root_idx in chunk_dirs:
        _remove_path(out_root_idx)


def _run_loc2ssst_parallel(params, ssst_control_file_tmp, iteration):
    """Compute the SSST correction grids of one iteration with loc2ssstCores
    parallel Loc2ssst instances.

    Every instance reads the SAME event locations (LSLOCFILES in the shared
    control file) but computes/writes grids only for ITS subset of stations
    (LSSTATIONS). Each station's correction depends only on that station's
    residuals, so the split is pure parallelisation.
    """
    # --- station list split round-robin into disjoint sets (sizes differ by
    # at most 1); empty sets (fewer stations than cores) are filtered out
    codes = _read_station_codes(params.stationsFile)
    n = params.loc2ssstCores
    sta_sets = ['='.join(codes[i::n]) for i in range(n)]
    sta_sets = [s for s in sta_sets if s]

    procs = []
    for index, sta_set in enumerate(sta_sets):
        n_sta = sta_set.count('=') + 1
        _print(params, f'Running: Loc2ssst instance {index} ({n_sta} stations)')
        control_file_idx = os.path.join(params.tmpDir, f'ssst_{iteration}.in_{index}')
        shutil.copy2(ssst_control_file_tmp, control_file_idx)
        with open(control_file_idx, 'a') as f:
            f.write(f'LSSTATIONS {sta_set}\n')
        log_path = f'{control_file_idx}.log'
        log_file = open(log_path, 'w')
        proc = subprocess.Popen([_exe('Loc2ssst'), control_file_idx],
                                cwd=params.tmpDir,
                                stdout=log_file, stderr=subprocess.STDOUT)
        procs.append((proc, log_path, log_file))

    try:
        _wait_all(params, [(p, lp) for p, lp, _ in procs], 'Loc2ssst')
    finally:
        for _, _, log_file in procs:
            log_file.close()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_grids(params, rebuild=False):
    """Build the initial P and S velocity/travel-time grids for the zone.

    Two full passes, one per phase: Vel2Grid + Grid2Time with the NLLoc
    control file as generated (VGTYPE P / GTFILES ... P), then again with a
    temporary S variant of it. Skipped when both P and S time grids already
    exist (resume support), unless rebuild=True.
    """
    p_grids = glob.glob(f'{params.initTimeRoot}.P.*.time.hdr')
    s_grids = glob.glob(f'{params.initTimeRoot}.S.*.time.hdr')
    if p_grids and s_grids and not rebuild:
        _print(params, f'Grids already present at {params.initTimeRoot} '
                       f'({len(p_grids)} P, {len(s_grids)} S) - skipping build')
        return

    os.makedirs(params.tmpDir, exist_ok=True)
    _check_disk(os.path.dirname(params.initTimeRoot) or '.')

    # P pass: the control file as generated
    _run_checked(params, [_exe('Vel2Grid'), params.fileRunNLL], 'Vel2Grid (P)',
                 os.path.join(params.tmpDir, 'vel2grid_P.log'))
    _run_checked(params, [_exe('Grid2Time'), params.fileRunNLL], 'Grid2Time (P)',
                 os.path.join(params.tmpDir, 'grid2time_P.log'))

    # S pass: temp variant with VGTYPE and the GTFILES phase switched to S
    s_control = os.path.join(params.tmpDir, 'run_S_variant.in')
    with open(params.fileRunNLL) as f:
        lines = f.readlines()
    with open(s_control, 'w') as f:
        for line in lines:
            stripped = line.lstrip()
            if stripped.startswith('VGTYPE'):
                f.write('VGTYPE S\n')
            elif stripped.startswith('GTFILES'):
                f.write(re.sub(r'\sP\s*$', ' S\n', line))
            else:
                f.write(line)

    _run_checked(params, [_exe('Vel2Grid'), s_control], 'Vel2Grid (S)',
                 os.path.join(params.tmpDir, 'vel2grid_S.log'))
    _run_checked(params, [_exe('Grid2Time'), s_control], 'Grid2Time (S)',
                 os.path.join(params.tmpDir, 'grid2time_S.log'))


def run_ssst(params, iteration_start=0, iteration_stop=None):
    """Run the iterative NLLoc + Loc2ssst SSST relocation workflow.

    Parameters
    ----------
    params : SSSTRunParams
    iteration_start : int, optional
        First iteration to run (default 0). Iteration N uses charDists[N];
        when > 0, the run resumes from the SSST grids of iteration
        ``iteration_start - 1`` (vpvsIter applies).
    iteration_stop : int, optional
        Last iteration to run, inclusive (default: len(charDists), the final
        NLLoc-only relocation). Use a smaller value for a partial campaign,
        then resume later with ``iteration_start = iteration_stop + 1``.

    Returns
    -------
    dict
        ``output``       — final (or last executed) summary .hyp path,
        ``csv``          — matching summary .csv path,
        ``log``          — path to the run journal,
        ``n_iterations`` — number of NLLoc iterations executed,
        ``next_iteration`` — resume point for a partial run, None when the
        final relocation was executed.
    """
    iteration_final = len(params.charDists)
    if iteration_stop is None:
        iteration_stop = iteration_final
    if not 0 <= iteration_start <= iteration_final:
        raise ValueError(f'iteration_start must be in 0..{iteration_final}.')
    if not iteration_start <= iteration_stop <= iteration_final:
        raise ValueError(
            f'iteration_stop must be in {iteration_start}..{iteration_final}.')

    # --- single source of truth for the min-phases threshold: the LSPHSTAT
    # NRdgsMin value of the Loc2ssst control file
    min_num_phases_loc = _read_min_num_phases(params.fileRunSSST)
    _print(params, f'Min number of phases (LSPHSTAT NRdgsMin): {min_num_phases_loc}')

    # --- initial state: iteration 0 locates with the original Grid2Time
    # grids and vpvs; resuming at N > 0 uses iteration N-1's corrected grids
    iteration = iteration_start
    iteration_max = iteration_final - 1
    if iteration == 0:
        time_root = params.initTimeRoot
        vpvs = params.vpvs
    else:
        time_root = os.path.join(
            params.outRoot, params.runName, params.ssstModelName,
            f'ssst_corr{iteration - 1}', params.catalog, params.model)
        vpvs = params.vpvsIter

    # --- fail-fast on missing grids at this run's actual starting point
    # (never built / wrong root on a fresh start, stale iteration_start on
    # a resume)
    if not glob.glob(f'{time_root}*.time.*'):
        raise FileNotFoundError(
            f'no travel-time grids match: {time_root}*.time.* ' +
            ('(run build_grids() first)' if iteration == 0 else
             '(check iteration_start against the last completed iteration)'))

    # --- archive the run configuration (both control files) so the run can
    # later be reproduced/understood; the journal joins them at the end
    out_run_dir = os.path.join(params.outRoot, params.runName)
    os.makedirs(out_run_dir, exist_ok=True)
    shutil.copy2(params.fileRunNLL, out_run_dir)
    shutil.copy2(params.fileRunSSST, out_run_dir)
    os.makedirs(params.tmpDir, exist_ok=True)
    os.makedirs(os.path.dirname(params.logFile), exist_ok=True)

    # --- run journal: fresh on a full run, appended on a resume
    if iteration_start == 0:
        with open(params.logFile, 'w') as f:
            f.write(f'{params.projectName}/{params.runName}/'
                    f'{params.ssstModelName} {params.catalog}\n')
    else:
        _log_append(params, f'Resumed at iteration {iteration_start}\n')

    # --- shared NLLoc control-file skeleton; every NLLoc instance of every
    # iteration starts from a copy of it
    skeleton_conf = _build_skeleton_conf(params)

    _print(params, f'Running: iterations {iteration}..{iteration_stop} '
                   f'(final={iteration_final}), '
                   f'CharDists={params.charDists[iteration:iteration_stop + 1]}')

    n_iterations = 0
    out_root = None
    ran_final = False

    # ======================== main iteration loop ==========================
    # each pass = one NLLoc relocation, followed (except on the final pass)
    # by one Loc2ssst correction computation feeding the next pass
    while iteration <= iteration_stop:
        _print(params, f'Start iteration {iteration}/{iteration_max} '
                       '=================================')
        _log_append(params, f'Start iteration {iteration}/{iteration_max} '
                            '=================================\n')

        # ===================== NLLoc step (every iteration) ================
        loc_obs = params.locObs
        save_nlloc_octree = params.saveNllocOctree
        save_fmamp = params.saveFmamp

        out_root = os.path.join(params.outRoot, params.runName,
                                params.ssstModelName,
                                f'loc_ssst_corr{iteration}', params.catalog)
        if iteration == iteration_final:
            # final iteration: relocate the FULL catalog (catalogFinal,
            # possibly a superset) with the last SSST grids, saving the
            # oct-tree grids + fmamp output for later processing
            out_root = os.path.join(params.outRoot, params.runName,
                                    params.ssstModelName,
                                    f'loc_ssst_corr{iteration}',
                                    params.catalogFinal)
            loc_obs = params.locObsFinal
            save_nlloc_octree = 'SAVE_NLLOC_OCTREE'
            save_fmamp = 'SAVE_FMAMP'

        os.makedirs(out_root, exist_ok=True)
        _log_append(params, f'    Running NLLoc: {out_root}\n')

        # per-iteration copy of the skeleton; _run_nlloc_parallel derives
        # the per-chunk control files from it
        control_file_tmp = os.path.join(params.tmpDir,
                                        f'{params.projectName}_nll_{iteration}.in')
        shutil.copy2(skeleton_conf, control_file_tmp)

        _run_nlloc_parallel(params, control_file_tmp, out_root, loc_obs,
                            time_root, min_num_phases_loc, vpvs,
                            save_nlloc_octree, save_fmamp, iteration)

        # from here on out_root is the FILE ROOT of this iteration's
        # locations (e.g. .../loc_ssst_corr2/GLOBAL_1/Pyrenees_1)
        out_root = os.path.join(out_root, params.projectName)
        n_iterations += 1

        # after the final NLLoc-only relocation there is no Loc2ssst step
        if iteration == iteration_final:
            ran_final = True
            break

        # ================= Loc2ssst step (SSST iterations only) =============
        char_dist = params.charDists[iteration]
        _print(params, f'Running Loc2ssst for iteration {iteration}/{iteration_max}  '
                       f'CharDist={char_dist} =================================')

        # output directory for this iteration's corrected travel-time grids
        ls_out_root = os.path.join(params.outRoot, params.runName,
                                   params.ssstModelName,
                                   f'ssst_corr{iteration}', params.catalog)
        os.makedirs(ls_out_root, exist_ok=True)

        # symlink ALL current travel-time grids into the new directory
        # first: Loc2ssst overwrites the links of the stations it corrects
        # with real files, while stations without corrections keep resolving
        # to their previous grid instead of disappearing
        for time_file in glob.glob(f'{time_root}*.time.*'):
            link = os.path.join(ls_out_root, os.path.basename(time_file))
            try:
                os.symlink(time_file, link)
            except FileExistsError:
                pass  # rerunning an iteration: links already in place

        # like out_root above, ls_out_root becomes a file root
        ls_out_root = os.path.join(ls_out_root, params.model)
        _log_append(params, f'    Running Loc2ssst  CharDist={char_dist} {ls_out_root}\n')

        # shared Loc2ssst control file of this iteration = fileRunSSST +
        # the iteration-specific statements:
        #   LSPARAMS   - smoothing distance (km) + weight floor
        #   LSMODE     - ANGLES_NO: no take-off angle grid regeneration
        #   LSOUT      - output root for correction/corrected grids
        #   LSLOCFILES - per-event .hyp locations of the NLLoc step above
        #   LOCFILES   - obs pattern + INPUT time-grid root corrected onto
        #   LOCMETH    - same location method settings as the NLLoc step
        # (_run_loc2ssst_parallel appends the per-instance LSSTATIONS line)
        ssst_control_file_tmp = os.path.join(
            params.tmpDir, f'{params.projectName}_ssst_{iteration}.in')
        shutil.copy2(params.fileRunSSST, ssst_control_file_tmp)
        with open(ssst_control_file_tmp, 'a') as f:
            f.write(f'LSPARAMS {char_dist} 0.0000001\n')
            f.write('LSMODE ANGLES_NO\n')
            f.write(f'LSOUT {ls_out_root}\n')
            f.write(f'LSLOCFILES {out_root}.*.*.grid0.loc.hyp\n')
            f.write(f'LOCFILES {loc_obs} NLLOC_OBS  {time_root}  {out_root}  0\n')
            f.write(f'LOCMETH EDT_OT_WT 9999.0 {min_num_phases_loc} -1 -1 {vpvs} 145 -1 1\n')

        _run_loc2ssst_parallel(params, ssst_control_file_tmp, iteration)

        # chain the iterations: the grids just written become the input
        # travel-time grids of the next iteration's NLLoc step
        time_root = ls_out_root

        _print(params, f'Finished iteration {iteration}/{iteration_max} '
                       '=================================')

        # once corrections are active, the (possibly different) Vp/Vs ratio
        # for corrected grids applies to all subsequent iterations
        iteration += 1
        vpvs = params.vpvsIter

    # ========================== end of the loop ============================
    # the run is complete only when the final NLLoc-only relocation actually
    # executed (a partial run stopping at the last SSST iteration still ends
    # with iteration == iteration_final without having run it)
    finished = ran_final
    next_iteration = None if finished else iteration

    if finished:
        ssst_model_dir = os.path.join(params.outRoot, params.runName,
                                      params.ssstModelName)
        # suggested manual cleanup once the run is validated: intermediate
        # grids (all but the last ssst_corr) and intermediate locations
        # (all but the final loc_ssst_corr) are deletable. NOT executed here
        # (the symlink chains make automatic deletion unsafe).
        k = iteration_final
        _log_append(params, '\n# cleanup commands (after validation only):\n')
        if k >= 2:
            _log_append(params, f'# rm -r {ssst_model_dir}/ssst_corr{{0..{k - 2}}}\n')
        _log_append(params, f'# rm -r {ssst_model_dir}/loc_ssst_corr{{0..{k - 1}}}\n')

        # record the definitive catalog in the ssst list
        with open(params.ssstListFile, 'a') as f:
            f.write(f'{out_root}.sum.grid0.loc.hyp\n')
    else:
        _log_append(params, f'Stopped after iteration {iteration_stop}; '
                            f'resume with iteration_start={next_iteration}\n')
        _print(params, f'Partial run done - resume with '
                       f'iteration_start={next_iteration}')

    # archive the journal next to the control-file copies
    shutil.copy2(params.logFile, out_run_dir)

    return {
        'output':         f'{out_root}.sum.grid0.loc.hyp',
        'csv':            f'{out_root}.sum.grid0.loc.csv',
        'log':            params.logFile,
        'n_iterations':   n_iterations,
        'next_iteration': next_iteration,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Iterative NLLoc + Loc2ssst SSST relocation workflow for one zone.')
    parser.add_argument('--run-name',       required=True)
    parser.add_argument('--project-name',   required=True,
                        help='Basename of the NLLoc output files (e.g. Pyrenees_1)')
    parser.add_argument('--model',          required=True,
                        help='Basename of the velocity/time grids (e.g. Pyrenees_1)')
    parser.add_argument('--catalog',        required=True, help='e.g. GLOBAL_1')
    parser.add_argument('--control-nll',    required=True, help='NLLoc control file')
    parser.add_argument('--control-ssst',   required=True, help='Loc2ssst control file')
    parser.add_argument('--loc-obs',        required=True,
                        help="Glob of per-event obs files (quote it: 'dir/*.nlloc_obs')")
    parser.add_argument('--stations',       required=True, help='GTSRCE station file')
    parser.add_argument('--init-time-root', required=True,
                        help='Initial travel-time grid file root')
    parser.add_argument('--out',            required=True, help='Output root directory')
    parser.add_argument('--tmp-dir',        required=True, help='Scratch directory')
    parser.add_argument('--log-file',       default=None,
                        help='Run journal path (default: <tmp-dir>/run_ssst.log)')
    parser.add_argument('--ssst-list',      default=None,
                        help='ssst.list path (default: <tmp-dir>/ssst.list)')
    parser.add_argument('--char-dists',     default='9999,50,15,5,1',
                        help='Comma-separated smoothing distances in km')
    parser.add_argument('--vpvs',           type=float, default=-9.99)
    parser.add_argument('--vpvs-iter',      type=float, default=-9.99)
    parser.add_argument('--nlloc-cores',    type=int, default=10)
    parser.add_argument('--loc2ssst-cores', type=int, default=10)
    parser.add_argument('--build-grids',    action='store_true',
                        help='Build the initial P and S grids before running')
    parser.add_argument('--rebuild-grids',  action='store_true',
                        help='Force the grid build even when grids exist')
    parser.add_argument('--iteration-start', type=int, default=0)
    parser.add_argument('--iteration-stop',  type=int, default=None,
                        help='Last iteration to run, inclusive (default: final relocation)')
    args = parser.parse_args()

    params = SSSTRunParams(
        runName      = args.run_name,
        projectName  = args.project_name,
        model        = args.model,
        catalog      = args.catalog,
        fileRunNLL   = args.control_nll,
        fileRunSSST  = args.control_ssst,
        locObs       = args.loc_obs,
        stationsFile = args.stations,
        initTimeRoot = args.init_time_root,
        outRoot      = args.out,
        tmpDir       = args.tmp_dir,
        logFile      = args.log_file or os.path.join(args.tmp_dir, 'run_ssst.log'),
        ssstListFile = args.ssst_list or os.path.join(args.tmp_dir, 'ssst.list'),
        charDists    = [float(v) if '.' in v else int(v)
                        for v in args.char_dists.split(',')],
        vpvs          = args.vpvs,
        vpvsIter      = args.vpvs_iter,
        nllocCores    = args.nlloc_cores,
        loc2ssstCores = args.loc2ssst_cores,
    )

    if args.build_grids or args.rebuild_grids:
        build_grids(params, rebuild=args.rebuild_grids)

    run_ssst(params, iteration_start=args.iteration_start,
             iteration_stop=args.iteration_stop)


if __name__ == '__main__':
    main()
