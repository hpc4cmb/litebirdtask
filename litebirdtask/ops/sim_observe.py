# Copyright (c) 2015-2023 by the parties listed in the AUTHORS file.
# All rights reserved.  Use of this source code is governed by
# a BSD-style license that can be found in the LICENSE file.

from datetime import datetime, timezone

import traitlets
import numpy as np

from astropy import units as u
from astropy.table import Column, QTable

import toast.qarray as qa
from toast.traits import trait_docs, Int, Unicode, Float, Bool, Instance, Quantity, Unit
from toast.utils import Environment, Logger
from toast.dist import distribute_discrete
from toast.timing import function_timer, Timer
from toast import Telescope, SpaceSite
from toast.ops import Operator, SimSatellite
from toast.ops.sim_satellite import satellite_scanning
from toast.ops.sim_hwp import simulate_hwp_response
from toast.observation import Observation, Session
from toast.observation import default_values as defaults
from toast.schedule import SatelliteSchedule

from ..instrument import load_imo


@trait_docs
class SimObserve(Operator):
    """Simulate satellite scanning.

    This operator loads a set of detectors from the IMO and observes with the specified
    schedule.

    """

    # Class traits (we also inherit the traits of the toast.ops.SimSatellite class)

    imo_file = Unicode(None, allow_none=True, help="The path to an IMO file.")

    select_telescope = Unicode(None, allow_none=True, help="Regular expression matching on telescope name")

    select_channel = Unicode(None, allow_none=True, help="Regular expression matching on channel name")
    
    select_wafer = Unicode(None, allow_none=True, help="Regular expression matching on wafer name")
    
    select_pixel = Unicode(None, allow_none=True, help="Regular expression matching on pixel name")
    
    select_detector = Unicode(None, allow_none=True, help="Regular expression matching on detector name")

    schedule = Instance(
        klass=SatelliteSchedule, allow_none=True, help="Instance of a SatelliteSchedule"
    )

    detset_key = Unicode(
        None,
        allow_none=True,
        help="If specified, use this column of the focalplane detector_data to group detectors",
    )

    times = Unicode(defaults.times, help="Observation shared key for timestamps")

    shared_flags = Unicode(
        defaults.shared_flags,
        allow_none=True,
        help="Observation shared key for common flags",
    )

    hwp_angle = Unicode(
        None, allow_none=True, help="Observation shared key for HWP angle"
    )

    boresight = Unicode(
        defaults.boresight_radec, help="Observation shared key for boresight"
    )

    position = Unicode(defaults.position, help="Observation shared key for position")

    velocity = Unicode(defaults.velocity, help="Observation shared key for velocity")

    det_data = Unicode(
        defaults.det_data,
        allow_none=True,
        help="Observation detdata key to initialize",
    )

    det_data_units = Unit(
        defaults.det_data_units, help="Output units if creating detector data"
    )

    det_flags = Unicode(
        defaults.det_flags,
        allow_none=True,
        help="Observation detdata key for flags to initialize",
    )
    
    wafer_obs = Bool(False, help="If True, split detectors into observations by wafer, not channel")

    distribute_time = Bool(
        False,
        help="Distribute observation data along the time axis rather than detector axis",
    )

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

    def _exec(self, data, detectors=None, **kwargs):
        zaxis = np.array([0, 0, 1], dtype=np.float64)

        if detectors is not None:
            msg = "Detectors should be selected through the IMO, "
            msg += "not passed as a function argument."
            raise RuntimeError(msg)
    
        # Load the IMO
        obs_tel = None
        if data.comm.world_rank == 0:
            obs_tel = load_imo(
                self.imo_file,
                telescope=self.select_telescope,
                channel=self.select_channel,
                wafer=self.select_wafer,
                pixel=self.select_pixel,
                detector=self.select_detector,
                wafer_obs=self.wafer_obs,
            )
        if data.comm.comm_world is not None:
            obs_tel = data.comm.comm_world.bcast(obs_tel, root=0)

        print(f"IMO = {obs_tel}", flush=True)
        
        # Data distribution in the detector direction
        det_ranks = data.comm.group_size
        if self.distribute_time:
            det_ranks = 1

        # Distribute the schedule based on the time covered by each scan.
        # The global start is the beginning of the first scan.

        mission_start = self.schedule.scans[0].start

        if len(self.schedule.scans) == 0:
            raise RuntimeError("Schedule has no scans!")
        
        scan_seconds = list()
        for scan in self.schedule.scans:
            scan_seconds.append(
                int((scan.stop - scan.start).total_seconds())
            )

        print(f"mission start = {mission_start}")
        print(f"schedule = {self.schedule.scans}")
        print(f"schedule seconds = {scan_seconds}", flush=True)

        # Distribute the observing sessions uniformly among groups.  We take each scan and
        # weight it by the duration.

        groupdist = distribute_discrete(scan_seconds, data.comm.ngroups)
        group_firstobs = groupdist[data.comm.group][0]
        group_numobs = groupdist[data.comm.group][1]
        
        for obindx in range(group_firstobs, group_firstobs + group_numobs):
            scan = self.schedule.scans[obindx]
            ses_start = scan.start.timestamp()
            ses_stop = scan.stop.timestamp()

            session = Session(
                f"{scan.name}_{int(ses_start):10d}",
                start=scan.start.astimezone(timezone.utc),
                end=scan.stop.astimezone(timezone.utc),
            )

            # For this observing session, loop over telescope objects and create
            # observations.

            for tel in obs_tel:
                focalplane = tel.focalplane
                rate = focalplane.sample_rate.to_value(u.Hz)
                incr = 1.0 / rate

                ffirst = rate * (scan.start - mission_start).total_seconds()
                first = int(ffirst)
                if ffirst - first > 1.0e-3 * incr:
                    first += 1
                start = first * incr + mission_start.timestamp()
                scan_samples = 1 + int(rate * (scan.stop.timestamp() - start))
                stop = (scan_samples - 1) * incr + start
                print(f"{scan_samples} : {start:0.15e} {stop:0.15e}", flush=True)

                # Get the detector sets
                detsets = None
                if self.detset_key is not None:
                    detsets = focalplane.detector_groups(self.detset_key)

                ob = Observation(
                    data.comm,
                    tel,
                    scan_samples,
                    name=f"{scan.name}_{int(ses_start)}_{tel.name}",
                    session=session,
                    detector_sets=detsets,
                    process_rows=det_ranks,
                )

                # Create shared objects for timestamps, common flags, position,
                # and velocity.
                ob.shared.create_column(
                    self.times,
                    shape=(ob.n_local_samples,),
                    dtype=np.float64,
                )
                ob.shared.create_column(
                    self.shared_flags,
                    shape=(ob.n_local_samples,),
                    dtype=np.uint8,
                )
                ob.shared.create_column(
                    self.position,
                    shape=(ob.n_local_samples, 3),
                    dtype=np.float64,
                )
                ob.shared.create_column(
                    self.velocity,
                    shape=(ob.n_local_samples, 3),
                    dtype=np.float64,
                )

                # Rank zero of each grid column creates the data

                stamps = None
                position = None
                velocity = None
                q_prec = None

                if ob.comm_col_rank == 0:
                    stamps = np.linspace(
                        start,
                        stop,
                        num=ob.n_local_samples,
                        endpoint=True,
                        dtype=np.float64,
                    )

                    print(f"stamps = {stamps}", flush=True)

                    # Get the motion of the site for these times.
                    position, velocity = tel.site.position_velocity(stamps)

                    # Get the quaternions for the precession axis.  For now, assume that
                    # it simply points away from the solar system barycenter

                    pos_norm = np.sqrt((position * position).sum(axis=1)).reshape(-1, 1)
                    pos_norm = 1.0 / pos_norm
                    prec_axis = pos_norm * position
                    q_prec = qa.from_vectors(
                        np.tile(zaxis, ob.n_local_samples).reshape(-1, 3), prec_axis
                    )

                ob.shared[self.times].set(stamps, offset=(0,), fromrank=0)
                ob.shared[self.position].set(position, offset=(0, 0), fromrank=0)
                ob.shared[self.velocity].set(velocity, offset=(0, 0), fromrank=0)

                # Create boresight pointing

                satellite_scanning(
                    ob,
                    self.boresight,
                    sample_offset=first,
                    q_prec=q_prec,
                    spin_period=scan.spin_period,
                    spin_angle=tel.spin_boresight_angle,
                    prec_period=scan.prec_period,
                    prec_angle=95 * u.degree - tel.spin_boresight_angle,
                )

                # Set HWP angle

                simulate_hwp_response(
                    ob,
                    ob_time_key=self.times,
                    ob_angle_key=self.hwp_angle,
                    ob_mueller_key=None,
                    hwp_start=start * u.second,
                    hwp_rpm=tel.hwp_rpm,
                )

                # Optionally initialize detector data

                dets = ob.select_local_detectors(detectors)

                if self.det_data is not None:
                    exists_data = ob.detdata.ensure(
                        self.det_data,
                        dtype=np.float64,
                        detectors=dets,
                        create_units=self.det_data_units,
                    )

                if self.det_flags is not None:
                    exists_flags = ob.detdata.ensure(
                        self.det_flags, dtype=np.uint8, detectors=dets
                    )

                data.obs.append(ob)

    def _finalize(self, data, **kwargs):
        return

    def _requires(self):
        return dict()

    def _provides(self):
        prov = {
            "shared": [
                self.times,
                self.shared_flags,
                self.boresight,
                self.hwp_angle,
                self.position,
                self.velocity,
            ]
        }
        if self.det_data is not None:
            prov["detdata"].append(self.det_data)
        if self.det_flags is not None:
            prov["detdata"].append(self.det_flags)
        return prov
