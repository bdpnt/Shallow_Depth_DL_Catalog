"""
keep_awake.py
============================
Keep macOS from idle-sleeping while a long relocation run is in progress.

Starts `caffeinate -i -w <pid>` in the background: it holds an idle-sleep
assertion and exits by itself as soon as the calling process ends, whether the
run finished, failed, or was killed — nothing has to release it.

-i blocks idle sleep only. Closing the lid still sleeps a laptop; the run is
then frozen, not lost, and resumes on wake.

Used by run_NLL.py and run_SSST.py (disabled there with --decaff).
"""

import os
import shutil
import subprocess


def keep_awake():
    """Hold an idle-sleep assertion until this process exits.

    Returns the caffeinate process, or None when caffeinate is unavailable
    (non-macOS), in which case the run proceeds without it.
    """
    if shutil.which('caffeinate') is None:
        print('caffeinate not found - the machine may sleep during the run')
        return None
    return subprocess.Popen(['caffeinate', '-i', '-w', str(os.getpid())])
