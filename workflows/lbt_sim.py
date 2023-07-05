#!/usr/bin/env python3
# Copyright (c) 2023-2023 by enities listed in the top-level AUTHORS file.
# Full license can be found in the top level LICENSE file.

"""
This script creates simulated LiteBIRD TOAST detector timestreams and then uses them to
test analysis techniques.
"""

import argparse
import os
import sys
import traceback

import numpy as np
from astropy import units as u

from litebirdtask import workflow

import toast


def parse_config(operators, templates, comm):
    """Parse command line arguments and load any config files.

    Return the final config, remaining args, and job size args.

    """
    # Argument parsing
    parser = argparse.ArgumentParser(description="LiteBIRD Simulation")

    # Arguments for this workflow
    parser.add_argument(
        "--out_dir",
        required=False,
        type=str,
        default=".",
        help="Output directory.",
    )

    parser.add_argument(
        "--debug",
        required=False,
        default=False,
        action="store_true",
        help="Make additional plots / checks for debugging",
    )

    # Build a config dictionary starting from the operator defaults, overriding with any
    # config files specified with the '--config' commandline option, followed by any
    # individually specified parameter overrides.

    config, args, jobargs = toast.parse_config(
        parser,
        operators=operators,
        templates=templates,
    )

    # Create our output directory
    if comm is None or comm.rank == 0:
        if not os.path.isdir(args.out_dir):
            os.makedirs(args.out_dir)

    # Log the config that was actually used at runtime.
    outlog = os.path.join(args.out_dir, "config_log.toml")
    toast.config.dump_toml(outlog, config, comm=comm)

    return config, args, jobargs


@toast.timing.function_timer
def main():
    env = toast.utils.Environment.get()
    log = toast.utils.Logger.get()
    gt = toast.timing.GlobalTimers.get()
    gt.start("lbt_sim (total)")

    # Get optional MPI parameters
    comm, procs, rank = toast.get_world()

    # The operators we want to configure from the command line or a parameter file.
    # We will use other operators, but these are the ones that the user can configure.
    # The "name" of each operator instance controls what the commandline and config
    # file options will be called.
    #
    # We can also set some default values here for the traits, including whether an
    # operator is disabled by default.

    operators = list()
    templates = list()

    # Simulated observing
    workflow.add_sim_observe_operators(operators)

    # Simulated data
    workflow.add_sim_operators(operators)

    # Default noise model
    workflow.add_simple_noise_operators(operators)

    # Pointing model
    workflow.add_pointing_operators(operators)

    # Noise model
    workflow.add_noise_operators(operators)

    # Mapmaking
    workflow.add_mapmaking_operators(operators, templates)

    # Parse options
    config, args, jobargs = parse_config(operators, templates, comm)

    # Instantiate our operators / templates that were configured from the
    # command line / files
    job = toast.create_from_config(config)

    # Simulate observing

    data = workflow.sim_observe(job, args, comm)

    # Simulate data

    workflow.run_simple_models(job, args, data)
    workflow.sim_data(job, args, comm)

    # Processing

    workflow.select_pointing(job, args, data)
    workflow.run_noise_estimation(job, args, data)
    workflow.run_mapmaking(job, args, data)

    # Collect optional timing information
    alltimers = toast.timing.gather_timers(comm=data.comm.comm_world)
    if data.comm.world_rank == 0:
        out = os.path.join(args.out_dir, "timing")
        toast.timing.dump(alltimers, out)


if __name__ == "__main__":
    world, procs, rank = toast.mpi.get_world()
    with toast.mpi.exception_guard(comm=world):
        main()

