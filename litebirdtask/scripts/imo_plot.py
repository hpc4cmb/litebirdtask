# Copyright (c) 2015-2023 by enities listed in the top-level AUTHORS file.
# Full license can be found in the top level LICENSE file.
"""Plot an instrument model.
"""

import argparse

from ..instrument import load_imo

from ..vis import plot_detectors


def main():
    parser = argparse.ArgumentParser(
        description="This program reads an IMO and plots the detectors",
        usage="lbt_imo_plot [options] (use --help for details)",
    )

    # parser.add_argument(
    #     "--hardware", required=True, default=None, help="Input hardware file"
    # )

    # parser.add_argument(
    #     "--out", required=False, default=None, help="Name of the output PDF file."
    # )

    # parser.add_argument(
    #     "--width",
    #     required=False,
    #     default=None,
    #     help="The width of the plot in degrees.",
    # )

    # parser.add_argument(
    #     "--height",
    #     required=False,
    #     default=None,
    #     help="The height of the plot in degrees.",
    # )

    # parser.add_argument(
    #     "--labels",
    #     required=False,
    #     default=False,
    #     action="store_true",
    #     help="Add pixel and polarization labels to the plot.",
    # )

    # args = parser.parse_args()

    # outfile = args.out
    # if outfile is None:
    #     fields = args.hardware.split(".")
    #     outfile = fields[0]

    # print("Loading hardware file {}...".format(args.hardware), flush=True)
    # hw = Hardware(args.hardware)

    # print("Generating detector plot...", flush=True)
    # plot_detectors(
    #     hw.data["detectors"],
    #     outfile,
    #     width=args.width,
    #     height=args.height,
    #     labels=args.labels,
    # )

    return
