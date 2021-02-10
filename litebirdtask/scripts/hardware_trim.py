# Copyright (c) 2015-2021 LiteBIRD Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Trim a hardware model to include only some detectors.
"""

import argparse

from ..hardware import Hardware


def main():
    parser = argparse.ArgumentParser(
        description="This program reads a hardware model from disk,\
            selects a subset of detectors, and writes the new model out.",
        usage="lb_hardware_trim [options] (use --help for details)",
    )

    parser.add_argument(
        "--hardware", required=True, default=None, help="Input hardware file"
    )

    parser.add_argument(
        "--out",
        required=False,
        default="trimmed",
        help="Name (without extensions) of the output hardware file",
    )

    parser.add_argument(
        "--plain",
        required=False,
        default=False,
        action="store_true",
        help="Write plain text (without gzip compression)",
    )

    parser.add_argument(
        "--overwrite",
        required=False,
        default=False,
        action="store_true",
        help="Overwrite any existing output file.",
    )

    parser.add_argument(
        "--telescopes",
        required=False,
        default=None,
        help="Select only detectors on these telescope (LFT, MFT, HFT)\
            .  This should be either a regex string or a comma-separated\
            list of names.",
    )

    parser.add_argument(
        "--match",
        required=False,
        default=None,
        nargs="*",
        help="Specify one or more detector match criteria.  Each match should\
            have the format '<property>:<regex or list>'.  The regex\
            expression should be valid to pass to the 're' module.  If \
            passing a list, this should be comma-separated.  For example, \
            --match 'band:.*040' 'wafer:L00' 'pol:A' ",
    )

    args = parser.parse_args()

    telescopes = args.telescopes
    if telescopes is not None:
        telescopes = telescopes.split(",")

    match = dict()
    if args.match is not None:
        for amat in args.match:
            fields = amat.split(":")
            if len(fields) != 2:
                raise ValueError("Invalid match string- must have a colon separator.")
            temp = fields[1].split(",")
            if len(temp) > 1:
                match[fields[0]] = temp
            else:
                match[fields[0]] = fields[1]

    print("Loading hardware from {}...".format(args.hardware), flush=True)
    hw = Hardware(args.hardware)

    print("Selecting detectors from:")
    if telescopes is not None:
        telstr = ""
        for tele in telescopes:
            telstr = "{}'{}', ".format(telstr, tele)
        print("  telescopes = {}".format(telstr.rstrip(", ")), flush=True)
    for k, v in match.items():
        print("  {} = r'{}'".format(k, v), flush=True)

    newhw = hw.select(telescopes=telescopes, match=match)

    if args.plain:
        outpath = "{}.toml".format(args.out)
        print("Dumping selected config to {}...".format(outpath))
        newhw.dump(outpath, overwrite=args.overwrite, compress=False)
    else:
        outpath = "{}.toml.gz".format(args.out)
        print("Dumping selected config to {}...".format(outpath))
        newhw.dump(outpath, overwrite=args.overwrite, compress=True)

    return
