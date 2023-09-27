# Copyright (c) 2023-2023 by enities listed in the top-level AUTHORS file.
# Full license can be found in the top level LICENSE file.
"""High-level workflow data processing.
"""

import os

import numpy as np
from astropy import units as u

import toast
from toast.observation import default_values as defaults
import toast.ops


def add_simple_noise_operators(operators):
    """Add basic nominal noise models"""
    operators.append(toast.ops.DefaultNoiseModel(name="default_model"))


def add_noise_operators(operators):
    operators.append(
        toast.ops.NoiseEstim(
            name="noise_estim",
            out_model="noise_est",
            lagmax=1024,
            nbin_psd=128,
            nsum=1,
            naverage=64,
        )
    )
    operators.append(
        toast.ops.FitNoiseModel(
            name="fit_estim",
            out_model="fit_noise",
        )
    )
    operators.append(
        toast.ops.FlagNoiseFit(
            name="flag_fit",
            noise_model="fit_noise",
        )
    )


def add_mapmaking_operators(operators, templates):
    templates.append(
        toast.templates.Offset(
            name="baselines",
            use_noise_prior=False,
            step_time=2.0 * u.second,
        )
    )
    operators.append(toast.ops.BinMap(name="binner"))
    operators.append(
        toast.ops.MapMaker(
            name="litebird",
            solve_rcond_threshold=1e-3,
            map_rcond_threshold=1e-3,
            iter_min=30,
            iter_max=200,
            write_hits=True,
            write_map=True,
            write_binmap=True,
            write_cov=True,
        )
    )


def run_simple_models(job, args, data):
    log = toast.utils.Logger.get()
    world_comm = data.comm.comm_world

    # Configured operators for this job
    job_ops = job.operators

    # Timer for reporting the progress
    timer = toast.timing.Timer()
    timer.start()

    # Simple noise models
    job_ops.default_model.apply(data)

    # Set the default noise model to the nominal one.
    job.map_noise_model = job_ops.default_model.noise_model

    timer.stop()


def run_noise_estimation(job, args, data):
    log = toast.utils.Logger.get()
    world_comm = data.comm.comm_world

    # Configured operators for this job
    job_ops = job.operators

    # Timer for reporting the progress
    timer = toast.timing.Timer()
    timer.start()

    # Estimate noise.  If the noise estimation is disabled at runtime, use the
    # default noise model.

    if job_ops.noise_estim.enabled:
        job_ops.noise_estim.det_data = defaults.det_data
        job_ops.noise_estim.apply(data)
        log.info_rank("Estimated noise in", comm=world_comm, timer=timer)

        # Create a fit to this noise model
        if job_ops.fit_estim.enabled:
            job_ops.fit_estim.noise_model = job_ops.noise_estim.out_model
            job_ops.fit_estim.apply(data)
            log.info_rank("Fit 1/f noise model in", comm=world_comm, timer=timer)
            job.map_noise_model = job_ops.fit_estim.out_model
            # Flag detector outliers
            job_ops.flag_fit.apply(data)
            log.info_rank("Flag noise model outliers in", comm=world_comm, timer=timer)
        else:
            job.map_noise_model = job_ops.noise_estim.out_model

        if args.debug:
            from toast.vis import plot_noise_estim

            # Plot the noise estimate fits
            nse_dir = os.path.join(args.out_dir, "noise_estim")
            if data.comm.world_rank == 0:
                os.makedirs(nse_dir)
            if data.comm.comm_world is not None:
                data.comm.comm_world.barrier()
            for ob in data.obs:
                est_model = ob[job_ops.noise_estim.out_model]
                if job_ops.fit_estim.enabled:
                    fit_model = ob[job_ops.fit_estim.out_model]
                    for det in ob.local_detectors:
                        fname = os.path.join(nse_dir, f"{ob.name}_{det}.pdf")
                        plot_noise_estim(
                            fname,
                            est_model.freq(det),
                            est_model.psd(det),
                            fit_freq=fit_model.freq(det),
                            fit_psd=fit_model.psd(det),
                            semilog=True,
                        )
                else:
                    for det in ob.local_detectors:
                        fname = os.path.join(nse_dir, f"{ob.name}_{det}.pdf")
                        plot_noise_estim(
                            fname,
                            est_model.freq(det),
                            est_model.psd(det),
                            semilog=True,
                        )
    timer.stop()


def run_mapmaking(job, args, data):
    log = toast.utils.Logger.get()
    world_comm = data.comm.comm_world

    # Configured operators for this job
    job_ops = job.operators

    # Configured templates
    job_tmpls = job.templates

    # Timer for reporting the progress
    timer = toast.timing.Timer()
    timer.start()

    # Set up the binning operator

    job_ops.binner.pixel_dist = "pixel_dist"
    job_ops.binner.pixel_pointing = job.pixels_solve
    job_ops.binner.stokes_weights = job.weights_solve
    job_ops.binner.noise_model = job.map_noise_model

    # Set up the template matrix

    job_tmpls.baselines.noise_model = job.map_noise_model
    tmatrix = toast.ops.TemplateMatrix(
        templates=[job_tmpls.baselines],
    )

    # Set up the mapmaker

    if args.debug:
        job_ops.litebird.keep_solver_products = True
    job_ops.litebird.binning = job_ops.binner
    job_ops.litebird.template_matrix = tmatrix
    job_ops.litebird.det_data = defaults.det_data
    job_ops.litebird.output_dir = args.out_dir
    job_ops.litebird.apply(data)
    log.info_rank("Finished map-making in", comm=world_comm, timer=timer)
