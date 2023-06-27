# Copyright (c) 2015-2023 by the parties listed in the AUTHORS file.
# All rights reserved.  Use of this source code is governed by
# a BSD-style license that can be found in the LICENSE file.

import os
from datetime import datetime

from astropy import units as u

import toast.ops
from toast.schedule_sim_satellite import create_satellite_schedule
from toast.instrument_sim import plot_focalplane
from toast.tests.mpi import MPITestCase
from toast.tests._helpers import create_comm, close_data

from .. import ops
from ..instrument import load_imo
from ._helpers import create_outdir, create_scanning


class SimObserveTest(MPITestCase):
    def setUp(self):
        fixture_name = os.path.splitext(os.path.basename(__file__))[0]
        self.outdir = create_outdir(self.comm, fixture_name)
        self.toastcomm = create_comm(self.comm)
        self.imo_path = None
        if "LB_IMO_FILE" in os.environ:
            self.imo_path = os.environ["LB_IMO_FILE"]

    def test_imo(self):
        if self.imo_path is None:
            print("'LB_IMO_FILE' not in environment, skipping tests")
            return
        
        imo = load_imo(
            self.imo_path,
            telescope="MFT",
            channel="M2-119",
        )
        print(imo)

        first_fp = imo[0].focalplane
        
        default_pixel_colors = {
            "LP1": (1.0, 0.0, 0.0, 0.2),
            "LP2": (0.4, 0.4, 1.0, 0.2),
            "LP3": (0.4, 0.4, 1.0, 0.2),
            "LP4": (0.4, 0.4, 1.0, 0.2),
            "MP1": (0.4, 1.0, 0.4, 0.2),
            "MP2": (0.4, 1.0, 0.4, 0.2),
            "HP1": (1.0, 0.4, 0.4, 0.2),
            "HP2": (1.0, 0.4, 0.4, 0.2),
            "HP3": (1.0, 0.4, 0.4, 0.2),
        }

        face_color = dict()
        pol_color = dict()
        for d in first_fp.detectors:
            detcolor = "black"
            if first_fp[d]["pol"] == "T":
                detcolor = (1.0, 0.0, 0.0, 1.0)
            if first_fp[d]["pol"] == "B":
                detcolor = (0.0, 0.0, 1.0, 1.0)
            pol_color[d] = detcolor
            face_color[d] = default_pixel_colors[first_fp[d]["pixtype"]]

        plot_focalplane(
            focalplane=imo[0].focalplane,
            width=25.0 * u.degree,
            height=25.0 * u.degree,
            face_color=face_color,
            pol_color=pol_color,
            show_labels=True,
            xieta=True,
            outfile=os.path.join(self.outdir, "imo_test.pdf")
        )

    def test_exec(self):
        if self.imo_path is None:
            print("'LB_IMO_FILE' not in environment, skipping tests")
            return
        
        data = create_scanning(
            self.comm,
            self.imo_path,
            tel="LFT",
            channel="L1-040",
            wafer="L00",
            session_per_group=1,
            session_time=10.0 * u.minute,
            prec_period=10 * u.minute,
            spin_period=1 * u.minute,
        )
        
        data.info()

        # # Scan fast enough to cover some sky in a short amount of time.  Reduce the
        # # angles to achieve a more compact hit map.
        # sim_sat = ops.SimSatellite(
        #     name="sim_sat",
        #     telescope=tele,
        #     schedule=sch,
        #     hwp_angle=defaults.hwp_angle,
        #     hwp_rpm=1.0,
        #     spin_angle=30.0 * u.degree,
        #     prec_angle=65.0 * u.degree,
        # )
        # sim_sat.apply(data)

        # # Plot some pointing
        # plotdetpointing = ops.PointingDetectorSimple(
        #     boresight=defaults.boresight_radec,
        #     quats="pquats",
        # )
        # plotdetpointing.apply(data)
        # if data.comm.world_rank == 0:
        #     n_debug = 10
        #     bquat = np.array(data.obs[0].shared[defaults.boresight_radec][0:n_debug, :])
        #     dquat = data.obs[0].detdata["pquats"][:, 0:n_debug, :]
        #     invalid = np.array(data.obs[0].shared[defaults.shared_flags][0:n_debug])
        #     invalid &= defaults.shared_mask_invalid
        #     valid = np.logical_not(invalid)
        #     outfile = os.path.join(self.outdir, "pointing.pdf")
        #     plot_projected_quats(
        #         outfile, qbore=bquat, qdet=dquat, valid=valid, scale=1.0
        #     )

        # # Expand pointing and make a hit map.
        # detpointing = ops.PointingDetectorSimple()
        # pixels = ops.PixelsHealpix(
        #     nest=True,
        #     create_dist="pixel_dist",
        #     detector_pointing=detpointing,
        # )
        # pixels.nside_submap = 2
        # pixels.nside = 8
        # pixels.apply(data)

        # build_hits = ops.BuildHitMap(pixel_dist="pixel_dist", pixels=pixels.pixels)
        # build_hits.apply(data)

        # # Plot the hits

        # hit_path = os.path.join(self.outdir, "hits.fits")
        # write_healpix_fits(data[build_hits.hits], hit_path, nest=pixels.nest)

        # if data.comm.world_rank == 0:
        #     set_matplotlib_backend()
        #     import matplotlib.pyplot as plt

        #     hits = hp.read_map(hit_path, field=None, nest=pixels.nest)
        #     outfile = os.path.join(self.outdir, "hits.png")
        #     hp.mollview(hits, xsize=1600, nest=True)
        #     plt.savefig(outfile)
        #     plt.close()
        
        close_data(data)


