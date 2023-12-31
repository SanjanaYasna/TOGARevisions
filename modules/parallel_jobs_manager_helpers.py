"""Helper function related to parallelisation."""
import time
import os
from modules.common import to_log

__author__ = "Bogdan M. Kirilenko"

ITER_DURATION = 60  # CESAR jobs check interval
NF_DIR_NAME = "nextflow_logs"


def monitor_jobs(jobs_managers, die_if_sc_1=False):
    """Monitor parallel jobs if many batches run simultaneously."""
    to_log(f"## Stated polling cluster jobs until they done")
    iter_num = 0
    while True:  # Run until all jobs are done (or crashed)
        all_done = True  # default val, re-define if something is not done
        for job_manager in jobs_managers:
            # check if each process is still running
            rc = job_manager.check_status()
            if rc is None:
                all_done = False
        if all_done:
            to_log("### CESAR jobs done ###")
            break
        else:
            to_log(f"Polling iteration {iter_num}; already waiting {ITER_DURATION * iter_num} seconds.")
            time.sleep(ITER_DURATION)
            iter_num += 1

    if any(jm.return_code != 0 for jm in jobs_managers) and die_if_sc_1 is True:
        # some para/nextflow job died: critical issue
        # if die_if_sc_1 is True: terminate the program
        err = "Error! Some para/nextflow processes died!"
        # TODO: think about the best error class
        raise AssertionError(err)


def get_nextflow_dir(proj_location, nf_dir_arg):
    """Define nextflow directory."""
    if nf_dir_arg is None:
        default_dir = os.path.join(proj_location, NF_DIR_NAME)
        os.mkdir(default_dir) if not os.path.isdir(default_dir) else None
        return default_dir
    else:
        os.mkdir(nf_dir_arg) if not os.path.isdir(nf_dir_arg) else None
        return nf_dir_arg
