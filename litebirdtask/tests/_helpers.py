# Copyright (c) 2015-2023 by the parties listed in the AUTHORS file.
# All rights reserved.  Use of this source code is governed by
# a BSD-style license that can be found in the LICENSE file.

import os
from datetime import datetime

from astropy import units as u

import toast
from toast.schedule_sim_satellite import create_satellite_schedule
from toast.tests._helpers import create_comm

from .. import ops


def create_outdir(mpicomm, subdir=None):
    """Create the top level output directory and per-test subdir.

    Args:
        mpicomm (MPI.Comm): the MPI communicator (or None).
        subdir (str): the sub directory for this test.

    Returns:
        str: full path to the test subdir if specified, else the top dir.

    """
    pwd = os.path.abspath(".")
    testdir = os.path.join(pwd, "litebirdtask_test_output")
    retdir = testdir
    if subdir is not None:
        retdir = os.path.join(testdir, subdir)
    if (mpicomm is None) or (mpicomm.rank == 0):
        if not os.path.isdir(testdir):
            os.mkdir(testdir)
        if not os.path.isdir(retdir):
            os.mkdir(retdir)
    if mpicomm is not None:
        mpicomm.barrier()
    return retdir


def create_scanning(
    mpicomm,
    imo_path,
    tel="LFT",
    channel="L1-040",
    wafer="L00",
    session_per_group=1,
    session_time=10.0 * u.minute,
    prec_period=10 * u.minute,
    spin_period=1 * u.minute,
):
    """Create a data object with a simulated scanning schedule.

    Use the specified MPI communicator to attempt to create 2 process groups.

    Args:
        mpicomm (MPI.Comm): the MPI communicator (or None).
        obs_per_group (int): the number of observations assigned to each group.
        obs_time (Quantity): the time length of one observation.

    Returns:
        toast.Data: the distributed data with named observations.

    """
    toastcomm = create_comm(mpicomm)
    data = toast.Data(toastcomm)

    # Create a schedule

    sch = create_satellite_schedule(
        prefix="test_",
        mission_start=datetime(2025, 2, 23),
        observation_time=session_time,
        gap_time=0 * u.minute,
        num_observations=(toastcomm.ngroups * session_per_group),
        prec_period=prec_period,
        spin_period=spin_period,
    )

    sim_obs = ops.SimObserve(
        imo_file=imo_path,
        select_telescope=tel,
        select_channel=channel,
        select_wafer=wafer,
        schedule=sch,
        detset_key="pixel",
    )
    sim_obs.apply(data)

    return data