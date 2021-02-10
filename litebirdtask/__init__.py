# Copyright (c) 2015-2021 LiteBIRD Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""LiteBIRD TOAST Analysis and Simulation Kit
"""

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# Namespace imports

from .hardware import Hardware
