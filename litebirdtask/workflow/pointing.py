# Copyright (c) 2021-2023 by enities listed in the top-level AUTHORS file.
# Full license can be found in the top level LICENSE file.
"""High-level workflow pointing operations.
"""

from astropy import units as u

import toast
from toast.observation import default_values as defaults


def add_pointing_operators(operators):
    operators.append(
        toast.ops.PointingDetectorSimple(
            name="det_pointing_radec", quats=defaults.quats
        )
    )
    operators.append(
        toast.ops.PixelsHealpix(
            name="pixels",
            nest=True,
            nside=1024,
            nside_submap=16,
        )
    )
    operators.append(
        toast.ops.PixelsHealpix(
            name="pixels_final",
            nest=True,
            nside=1024,
            nside_submap=16,
            enabled=False,
        )
    )
    operators.append(
        toast.ops.StokesWeights(
            name="weights_radec", weights="weights_radec", mode="IQU"
        )
    )


def select_pointing(job, args, data):
    """Select the pixelization scheme for both the solver and final binning."""
    log = toast.utils.Logger.get()

    job_ops = job.operators

    job_ops.det_pointing_radec.boresight = defaults.boresight_radec

    job_ops.pixels.detector_pointing = job_ops.det_pointing_radec

    job_ops.weights_radec.detector_pointing = job_ops.det_pointing_radec
    job_ops.weights_radec.hwp_angle = defaults.hwp_angle

    # Select Pixelization and weights for solve and final binning

    job.pixels_solve = job_ops.pixels
    job.weights_solve = job_ops.weights_radec

    if job_ops.pixels_final.enabled:
        # The user enabled a different pixelization for the final map
        job.pixels_final = job_ops.pixels_final
    else:
        job.pixels_final = job.pixels_solve
    job.weights_final = job.weights_solve

    log.info_rank(
        f"Template solve using pixelization: {job.pixels_solve.name}",
        comm=data.comm.comm_world,
    )
    log.info_rank(
        f"Template solve using weights: {job.weights_solve.name}",
        comm=data.comm.comm_world,
    )
    log.info_rank(
        f"Final binning using pixelization: {job.pixels_final.name}",
        comm=data.comm.comm_world,
    )
    log.info_rank(
        f"Final binning using weights: {job.weights_final.name}",
        comm=data.comm.comm_world,
    )
