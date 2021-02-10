# Copyright (c) 2015-2021 LiteBIRD Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Plot a hardware model.
"""

import argparse

from ..hardware import Hardware

from ..vis import summary_text


def main():
    parser = argparse.ArgumentParser(
        description="This program reads a hardware model and prints some\
            summary text to the terminal.",
        usage="lb_hardware_info [hardware file, [hardware_file]] ...",
    )

    parser.add_argument("hardware", type=str, nargs="+", help="Input hardware file")

    args = parser.parse_args()

    for file in args.hardware:
        print("Loading hardware file {}...".format(file), flush=True)
        hw = Hardware(file)
        summary_text(hw)

    return
