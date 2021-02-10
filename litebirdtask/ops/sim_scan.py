# Copyright (c) 2015-2021 LiteBIRD Collaboration.
# Full license can be found in the top level "LICENSE" file.

import traitlets

import numpy as np

from astropy import units as u

import toast.qarray as qa

from toast.traits import trait_docs, Int, Unicode, Float, Bool, Instance, Quantity

from toast.utils import Environment, name_UID, Logger, rate_from_times

from toast.dist import distribute_discrete

from ..timing import function_timer, Timer

from ..intervals import Interval

from toast.schedule import SatelliteSchedule

from toast import Observation, Telescope, Focalplane, SpaceSite

from toast.ops import Operator

from toast.ops import simulate_hwp_response

from toast.ops.sim_satellite import satellite_scanning, SimSatellite

import litebirdtask as lbt


@trait_docs
class SimScan(Operator):
    """Simulate satellite scanning.

    This uses an observing schedule to simulate the satellite motion.  The specified
    hardware object should already be trimmed to include just the detectors to be
    simulated.

    """

    # Class traits

    API = Int(0, help="Internal interface version for this operator")

    hardware = Instance(klass=lbt.Hardware, allow_none=True, help="The hardware model")

    schedule = Instance(
        klass=SatelliteSchedule, allow_none=True, help="Instance of a SatelliteSchedule"
    )

    spin_angle = Quantity(
        30.0 * u.degree, help="The opening angle of the boresight from the spin axis"
    )

    prec_angle = Quantity(
        65.0 * u.degree,
        help="The opening angle of the spin axis from the precession axis",
    )

    hwp_rpm = Float(None, allow_none=True, help="The rate (in RPM) of the HWP rotation")

    hwp_step = Quantity(
        None, allow_none=True, help="For stepped HWP, the angle of each step"
    )

    hwp_step_time = Quantity(
        None, allow_none=True, help="For stepped HWP, the time between steps"
    )

    distribute_time = Bool(
        False,
        help="Distribute observation data along the time axis rather than detector axis",
    )

    times = Unicode("times", help="Observation shared key for timestamps")

    shared_flags = Unicode("flags", help="Observation shared key for common flags")

    hwp_angle = Unicode("hwp_angle", help="Observation shared key for HWP angle")

    boresight = Unicode("boresight_radec", help="Observation shared key for boresight")

    position = Unicode("position", help="Observation shared key for position")

    velocity = Unicode("velocity", help="Observation shared key for velocity")

    @traitlets.validate("hardware")
    def _check_hardware(self, proposal):
        hw = proposal["value"]
        if hw is not None:
            if not isinstance(hw, lbt.Hardware):
                raise traitlets.TraitError(
                    "hardware must be a litebirdtask.Hardware instance"
                )
        return hw

    @traitlets.validate("schedule")
    def _check_schedule(self, proposal):
        sch = proposal["value"]
        if sch is not None:
            if not isinstance(sch, SatelliteSchedule):
                raise traitlets.TraitError(
                    "schedule must be an instance of a SatelliteSchedule"
                )
        return sch

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # For now, this class is almost identical to the one in TOAST.  We will just
        # use that class until something about LiteBIRD scanning diverges from that
        # of a generic satellite.  The corresponding instance of the TOAST class will
        # be instantiated on the first call to exec().
        self._sat_sim = None

    def _exec(self, data, detectors=None, **kwargs):
        zaxis = np.array([0, 0, 1], dtype=np.float64)
        log = Logger.get()
        comm = data.comm

        if self.hardware is None:
            raise RuntimeError(
                "The hardware attribute must be set before calling exec()"
            )
        if self.schedule is None:
            raise RuntimeError(
                "The schedule attribute must be set before calling exec()"
            )

        if self._sat_sim is None:
            # This is our first call to exec().  Create a toast Focalplane and Telescope
            # and then instantiate the generic toast SimSatellite class.

            site = SpaceSite("L2")

            focalplane = Focalplane()

            telescope = Telescope("LiteBIRD", focalplane=focalplane, site=site)

            self._sat_sim = SimSatellite(
                telescope=telescope,
                schedule=self.schedule,
                spin_angle=self.spin_angle,
                prec_angle=self.prec_angle,
                hwp_rpm=self.hwp_rpm,
                hwp_step=self.hwp_step,
                hwp_step_time=self.hwp_step_time,
                distribute_time=self.distribute_time,
                times=self.times,
                shared_flags=self.shared_flags,
                hwp_angle=self.hwp_angle,
                boresight=self.boresight,
                position=self.position,
                velocity=self.velocity,
            )

        # Run the generic satellite simulation
        self._sat_sim._exec(data, detectors=detectors, **kwargs)

        return

    def _finalize(self, data, **kwargs):
        if self._sat_sim is not None:
            # We have run exec() at least once...
            self._sat_sim._finalize(data, **kwargs)
        return

    def _requires(self):
        return dict()

    def _provides(self):
        return {
            "shared": [
                self.times,
                self.shared_flags,
                self.boresight,
                self.hwp_angle,
                self.position,
                self.velocity,
            ]
        }

    def _accelerators(self):
        return list()
