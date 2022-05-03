# Copyright (c) 2015-2021 LiteBIRD Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Export a pre-trimmed hardware model into a toast focalplane.
"""

import os

import argparse

from ..hardware import Hardware
from toast.io import H5File

def main():
    parser = argparse.ArgumentParser(
        description="This program reads a hardware model from disk \
            and writes out a TOAST focalplane file.",
        usage="lbt_export_focalplane [options] (use --help for details)",
    )

    parser.add_argument(
        "--hardware", required=True, default=None, help="Input hardware file"
    )

    parser.add_argument(
        "--out",
        required=False,
        default="lb_focalplane.h5",
        help="Name of the output HDF5 file",
    )

    parser.add_argument(
        "--overwrite",
        required=False,
        default=False,
        action="store_true",
        help="Overwrite any existing output file.",
    )

    args = parser.parse_args()

    print("Loading hardware from {}...".format(args.hardware), flush=True)
    hw = Hardware(args.hardware)

    print("Exporting Focalplane", flush=True)
    fp = hw.focalplane()

    if os.path.isfile(args.out):
        if args.overwrite:
            print(f"Overwriting output file {args.out}")
            os.remove(args.out)
        else:
            raise RuntimeError(f"File {args.out} exists and overwrite is False")
    else:
        print(f"Creating output file {args.out}")

    with H5File(args.out, "w",  force_serial=True) as f:
            fp.save_hdf5(f.handle)
    

