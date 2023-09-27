# Copyright (c) 2023-2023 by enities listed in the top-level AUTHORS file.
# Full license can be found in the top level LICENSE file.
"""High-level workflow data I/O.
"""

import os
from ast import literal_eval as make_tuple

import numpy as np
from astropy import units as u

import toast
from toast.observation import default_values as defaults
from toast.utils import Logger


def add_data_args(parser):
    """Add commandline args for data loading.

    Since there is no real data yet, this currently just supports loading
    toast Observation dumps.
    
    """
    parser.add_argument(
        "--obs_hdf5",
        required=False,
        action="extend",
        nargs="*",
        help="Path to a TOAST hdf5 observation dump (can use multiple times)",
    )


def load_data(job, args, mpi_comm):
    """Load one or more observations from disk.
    """
    log = Logger.get()

    # Configured operators for this job
    job_ops = job.operators

    # Create the toast communicator.  We will use one group
    # containing all processes.
    gsize = 1
    if mpi_comm is not None:
        gsize = mpi_comm.size
    toast_comm = toast.Comm(world=mpi_comm, groupsize=gsize)

    # The data container
    data = toast.Data(comm=toast_comm)

    for hobs in list(args.obs_hdf5):
        log.info_rank(f"Starting load of HDF5 data {hobs}", comm=data.comm.comm_group)
        ob = toast.io.load_hdf5(
            hobs,
            toast_comm,
            process_rows=toast_comm.group_size,
            meta=None,
            detdata=None,
            shared=None,
            intervals=None,
            force_serial=True,
        )
        data.obs.append(ob)

    return data

    