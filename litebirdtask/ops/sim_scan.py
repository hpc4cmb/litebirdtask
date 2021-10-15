# Copyright (c) 2015-2021 LiteBIRD Collaboration.
# Full license can be found in the top level "LICENSE" file.

import traitlets

import numpy as np

from astropy import units as u
from astropy.table import Column, QTable

from toast.traits import trait_docs, Int, Unicode, Float, Bool, Instance, Quantity

from toast.utils import Environment, Logger

from toast.timing import function_timer, Timer

from toast import Telescope, SpaceSite

from toast.ops.sim_satellite import SimSatellite

from ..hardware import Hardware


@trait_docs
class SimScan(SimSatellite):
    """Simulate satellite scanning.

    This uses an observing schedule to simulate the satellite motion.  The specified
    hardware object should already be trimmed to include just the detectors to be
    simulated.

    """

    # Class traits (we also inherit the traits of the toast.ops.SimSatellite class)

    hardware = Instance(klass=Hardware, allow_none=True, help="The hardware model")

    @traitlets.validate("hardware")
    def _check_hardware(self, proposal):
        hw = proposal["value"]
        if hw is not None:
            if not isinstance(hw, Hardware):
                raise traitlets.TraitError(
                    "hardware must be a litebirdtask.Hardware instance"
                )
            # Set up the telescope trait in the base class.
            site = SpaceSite("L2")
            self.telescope = Telescope("LiteBIRD", focalplane=hw.focalplane(), site=site)
        return hw

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _exec(self, data, detectors=None, **kwargs):
        super()._exec(data, detectors=detectors, **kwargs)

    def _finalize(self, data, **kwargs):
        return super()._finalize(data, **kwargs)

    def _requires(self):
        return super()._requires()

    def _provides(self):
        return super()._provides()

    def _accelerators(self):
        return list()
