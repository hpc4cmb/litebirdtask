# Copyright (c) 2021-2023 by enities listed in the top-level AUTHORS file.
# Full license can be found in the top level LICENSE file.
"""High-level workflow simulation functions.
"""

import os
import numpy as np

import toast
import toast.ops
from toast.observation import default_values as defaults

from .. import ops as lbtops


def add_sim_observe_operators(operators):
    """Add the litebird operator for simulated observing.
    """
    operators.append(
        lbtops.SimObserve(
            name="sim_observe",
            detset_key="pixel",
        )
    )


def add_sim_operators(operators):
    """Our standard set of simulation operators to configure from a job.
    """
    operators.extend(
        [
            toast.ops.SimDipole(
                name="sim_dipole",
                mode="total",
            ),
            toast.ops.ScanHealpixMap(name="scan_map"),
            toast.ops.SimNoise(name="sim_noise"),
            toast.ops.SaveHDF5(name="save_hdf5"),
        ]
    )


def sim_observe(job, args, data):
    """Run simulated observing.
    """
    ops = job.operators
    if ops.sim_observe.imo_file is None:
        msg = f"You must set the {ops.sim_observe.name}.imo_file trait before running"
        raise RuntimeError(msg)
    ops.sim_observe.apply(data)


def sim_data(job, args, data):
    """Perform nominal data simulation with configured operators."""
    log = toast.utils.Logger.get()
    ops = job.operators
    world_comm = data.comm.comm_world

    # Timer for reporting the progress
    timer = toast.timing.Timer()
    timer.start()

    # Set the detector data to zero if it is not already
    toast.ops.Reset(detdata=[defaults.det_data]).apply(data)

    # Load simulated sky signal
    if ops.scan_map.enabled and ops.scan_map.file is not None:
        log.info_rank("Begin simulated sky signal", comm=world_comm)
        ops.scan_map.apply(data)
        log.info_rank("Finished simulated sky signal in", comm=world_comm, timer=timer)

    # Simulated Dipole
    if ops.sim_dipole.enabled:
        log.info_rank("Begin dipole simulation", comm=world_comm)
        ops.sim_dipole.det_data = defaults.det_data
        ops.sim_dipole.apply(data)
        log.info_rank("Finished dipole sim in", comm=world_comm, timer=timer)

    # Accumulate instrumental noise
    if ops.sim_noise.enabled:
        log.info_rank("Begin instrument noise simulation", comm=world_comm)
        ops.sim_noise.noise_model = ops.default_model.out_model
        ops.sim_noise.det_data = defaults.det_data
        ops.sim_noise.apply(data)
        log.info_rank("Finished instrument noise sim in", comm=world_comm, timer=timer)

    # Optionally write out the data
    if ops.save_hdf5.enabled:
        log.info_rank("Begin save HDF5 data", comm=world_comm)
        if ops.save_hdf5.volume is None:
            ops.save_hdf5.volume = os.path.join(args.out_dir, "data")
        ops.save_hdf5.apply(data)
        log.info_rank("Finished save HDF5 data in", comm=world_comm, timer=timer)
