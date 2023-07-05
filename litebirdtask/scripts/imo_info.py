# Copyright (c) 2015-2023 by enities listed in the top-level AUTHORS file.
# Full license can be found in the top level LICENSE file.
"""Print information about an instrument model.
"""

import argparse

from ..instrument import load_imo

from ..vis import summary_text


def main():
    parser = argparse.ArgumentParser(
        description="This program reads an IMO and prints some\
            summary text to the terminal.",
        usage="lbt_imo_info [options] <IMO JSON file>",
    )

    parser.add_argument("imo", type=str, help="IMO file")

    parser.add_argument(
        "telescope", 
        type=str,
        default=None,
        help="Telescope selection as a regular expression"
    )

    parser.add_argument(
        "channel", 
        type=str,
        default=None,
        help="Channel selection as a regular expression"
    )

    parser.add_argument(
        "wafer", 
        type=str,
        default=None,
        help="Wafer selection as a regular expression"
    )

    parser.add_argument(
        "detector", 
        type=str,
        default=None,
        help="detector selection as a regular expression"
    )

    args = parser.parse_args()

    scan_props, tele_list = load_imo(
        args.imo,
        telescope=args.telescope,
        channel=args.channel,
        wafer=args.wafer,
        pixel=None,
        detector=args.detector,
    )

    summary_text(scan_props, tele_list)

    if len(tele_list) == 1:
        fp = tele_list[0].focalplane
        det_props = fp.detector_data
        if len(det_props) == 1:
            for col in det_props.colnames:
                print(f"  {col} = {det_props[col][0]}")

    return
