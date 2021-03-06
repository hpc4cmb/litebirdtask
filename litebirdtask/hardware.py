# Copyright (c) 2015-2021 LiteBIRD collaboration.
# Full license can be found in the top level "LICENSE" file.
"""Hardware configuration utilities.
"""

import os
import re
import copy

from collections import OrderedDict

import gzip

import toml


class Hardware(object):
    """Class representing a specific hardware configuration.

    The data is stored in a dictionary, and can be loaded / dumped to disk
    as well as trimmed to include only a subset of detectors.

    Args:
        path (str, optional): If specified, configuration is loaded from this
            file during construction.

    """

    def __init__(self, path=None):
        self.data = OrderedDict()
        if path is not None:
            self.load(path)

    def dump(self, path, overwrite=False, compress=False):
        """Write hardware config to a TOML file.

        Dump data to a TOML format file, optionally compressing the contents
        with gzip and optionally overwriting the file.

        Args:
            path (str): The file to write.
            overwrite (bool): If True, overwrite the file if it exists.
                If False, then existing files will cause an exception.
            compress (bool): If True, compress the data with gzip on write.

        Returns:
            None

        """
        if os.path.exists(path):
            if overwrite:
                os.remove(path)
            else:
                raise RuntimeError(
                    "Dump path {} already exists.  Use " "overwrite option".format(path)
                )
        if compress:
            with gzip.open(path, "wb") as f:
                dstr = toml.dumps(self.data)
                f.write(dstr.encode())
        else:
            with open(path, "w") as f:
                dstr = toml.dumps(self.data)
                f.write(dstr)
        return

    def load(self, path):
        """Read data from a TOML file.

        The file can either be regular text or a gzipped version of a TOML
        file.

        Args:
            path (str): The file to read.

        Returns:
            None

        """
        dstr = None
        try:
            with gzip.open(path, "rb") as f:
                dstr = f.read()
                self.data = toml.loads(dstr.decode())
        except OSError:
            with open(path, "r") as f:
                dstr = f.read()
                self.data = toml.loads(dstr)
        return

    def wafer_map(self):
        """Construct wafer mapping to other auxilliary data.

        Given the current data state, build dictionaries to go from wafers to all other
        non-detector info:  telescopes, bands, etc.  This is a convenient mapping when
        pruning the hardware information or doing other kinds of lookups.

        Returns:
            (dict): Nested dictionaries from wafers to other properties.

        """
        result = OrderedDict()

        wafer_to_tele = dict()
        for tb, props in self.data["telescopes"].items():
            for wf in props["wafers"]:
                wafer_to_tele[wf] = tb

        result["pixels"] = dict()
        result["bands"] = dict()
        for wafer_name, wafer_props in self.data["wafers"].items():
            result["pixels"][wafer_name] = wafer_props["pixels"]
            result["bands"][wafer_name] = list()
            for pix in result["pixels"][wafer_name]:
                result["bands"][wafer_name].extend(self.data["pixels"][pix]["bands"])

        result["telescopes"] = {
            x: wafer_to_tele[x] for x in list(self.data["wafers"].keys())
        }
        return result

    def select(self, telescopes=None, match=dict()):
        """Select a subset of detectors.

        Select detectors whose properties match some criteria.  A new Hardware
        object is created and returned.  If a matching expression is not
        specified for a given property name, then this is equivalent to
        selecting all values of that property.

        Before selecting on detector properties, any telescope filtering
        criteria are first applied.

        Each key of the "match" dictionary should be the name of a detector
        property to be considered for selection (e.g. band, wafer, pol, pixel).
        The value is a matching expression which can be:

            - A list of explicit values to match.
            - A string containing a regex expression to apply.

        Example:
            Imagine you wanted to select all 40GHz detectors on wafers L00 and
            L01 which have "A" polarization and are located in pixels 0-9
            (recall the "." matches a single character)::

                new = hw.select(match={"wafer": ["L00", "L01"],
                                "band": ".*040",
                                "pol": "A",
                                "pixel": "00."})

        Args:
            telescopes (str): A regex string to apply to telescope names or a
                list of explicit names.
            match (dict): The dictionary of property names and their matching
                expressions.

        Returns:
            (Hardware): A new Hardware instance with the selected detectors.

        """
        # First parse any telescope options into a list of wafers
        wselect = list()
        if telescopes is not None:
            for tele in telescopes:
                wselect.extend(self.data["telescopes"][tele]["wafers"])

        dets = self.data["detectors"]

        # Build regex matches for each property
        reg = dict()
        if "wafer" in match:
            # Handle wafer case separately, since we need to merge any
            # match with our telescope selection of wafers above.
            k = "wafer"
            v = match[k]
            if wselect is None:
                # Just the regular behavior
                if isinstance(v, list):
                    reg[k] = re.compile(r"(" + "|".join(v) + r")")
                else:
                    reg[k] = re.compile(v)
            else:
                # Merge our selection
                wall = list(wselect)
                if isinstance(v, list):
                    wall.extend(v)
                else:
                    wall.append(v)
                reg[k] = re.compile(r"(" + "|".join(wall) + r")")
        elif wselect is not None:
            # No pattern in the match dictionary, just our list from the
            # telescope selection.
            reg["wafer"] = re.compile(r"(" + "|".join(wselect) + r")")

        for k, v in match.items():
            if k == "wafer":
                # Already handled above
                continue
            else:
                if isinstance(v, list):
                    reg[k] = re.compile(r"(" + "|".join(v) + r")")
                else:
                    reg[k] = re.compile(v)

        # Go through all detectors selecting things that match all fields
        newwafers = set()
        newdets = OrderedDict()
        for d, props in dets.items():
            keep = True
            for k, v in reg.items():
                if k in props:
                    test = v.match(props[k])
                    if test is None:
                        keep = False
                        break
            if keep:
                newwafers.add(props["wafer"])
                newdets[d] = copy.deepcopy(props)

        # Now compute the reduced set of auxilliary data needed for these
        # detectors.
        wafermap = self.wafer_map()

        # Copy this data
        hw = Hardware()
        hw.data = OrderedDict()
        hw.data["software"] = copy.deepcopy(self.data["software"])
        for k, v in wafermap.items():
            hw.data[k] = OrderedDict()
            tocopy = set()
            for wf in newwafers:
                if isinstance(v[wf], list):
                    for iv in v[wf]:
                        tocopy.add(iv)
                else:
                    tocopy.add(v[wf])
            for elem in tocopy:
                hw.data[k][elem] = copy.deepcopy(self.data[k][elem])

        # Copy over the wafer data
        hw.data["wafers"] = OrderedDict()
        for wf in newwafers:
            hw.data["wafers"][wf] = copy.deepcopy(self.data["wafers"][wf])

        # And the detectors...
        hw.data["detectors"] = newdets

        return hw
