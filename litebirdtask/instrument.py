# Copyright (c) 2015-2023 by the parties listed in the AUTHORS file.
# All rights reserved.  Use of this source code is governed by
# a BSD-style license that can be found in the LICENSE file.
"""Instrument Classes
"""

import json
import os
import sys
from datetime import datetime, timezone, timedelta
import re

import numpy as np

import astropy.units as u
from astropy.table import QTable, Column
from toast.instrument import Focalplane, SpaceSite, Telescope
from toast.schedule import SatelliteSchedule, SatelliteScan
from toast.utils import Logger


class LitebirdSite(SpaceSite):
    def __init__(
        self,
        name="LiteBIRD",
        **kwargs,
    ):
        super().__init__(
            name,
            **kwargs,
        )

    def _position(self, times):
        # Eventually we could insert simulated or real orbital dynamics
        return super()._position(times)

    def _velocity(self, times):
        # Eventually we could insert simulated or real orbital dynamics
        return super()._velocity(times)


def load_imo(
    path,
    telescope=None,
    channel=None,
    wafer=None,
    pixel=None,
    detector=None,
    wafer_obs=False,
):
    """Load a subset of an Instrument Model.

    This function loads an IMO file and selects a subset of detectors.  The resulting
    detectors are grouped either by channel or wafer, and placed into a Telescope
    instance suitable for creating Observations.

    Args:
        path (str):  The path to the IMO JSON file.
        telescope (str):  Only use detectors from this telescope.
        channel (str):  Only use this channel (frequency).  The string is interpreted
            as a regular expression.
        wafer (str):  Only use this wafer.  The string is interpreted as a regular
            expression.
        pixel (str):  Only use this pixel.  The string is interpreted as a regular
            expression.
        detector (str):  Only use this detector.  The string is interpreted as a regular
            expression.
        wafer_obs (bool):  If True, split detectors into observations according to
            wafer.  Default is to split by telescope-channel.

    Returns:
        (tuple):  A tuple containing the IMO scanning parameters and a list of 
            Telescope objects, each suitable for creating an Observation.

    """
    if telescope is None:
        telescope = ".*"
    if channel is None:
        channel = ".*"
    if wafer is None:
        wafer = ".*"
    if pixel is None:
        pixel = ".*"
    if detector is None:
        detector = ".*"
    tel_pat = re.compile(telescope)
    chan_pat = re.compile(channel)
    wafer_pat = re.compile(wafer)
    pixel_pat = re.compile(pixel)
    det_pat = re.compile(detector)

    # Load the full IMO
    with open(path, "r") as f:
        imo = json.load(f)

    # Extract all info from the "data_files".
    tel_keys = [
        "channel_names",
        "hwp_rpm",
        "sampling_rate_hz",
        "spin_boresight_angle_deg",
        "boresight_rotangle_deg",
        "spin_rotangle_deg",
    ]
    chan_keys = [
        "sampling_rate_hz",
        "net_detector_ukrts",
        "net_channel_ukrts",
        "fknee_mhz",
        "pol_sensitivity_detector_ukarcmin",
        "pol_sensitivity_channel_ukarcmin",
        "bandcenter_ghz",
        "bandwidth_ghz",
        "fwhm_arcmin",
        "fmin_hz",
        "alpha",
        "psd_thermalbath_corr_ukcmb^2",
        "psd_thermalbath_uncorr_ukcmb^2",
    ]
    scan_keys = [
        "spin_sun_angle_deg",
        "precession_period_min",
        "spin_rate_rpm",
        "mission_duration_year",
        "observation_duty_cycle",
        "cosmic_ray_loss",
        "margin",
        "detector_yield",
    ]
    # Do one pass to extract all telescope and channel info
    tel_data = dict()
    scan_data = None
    chan_data = dict()
    chan_to_tel = dict()
    for iobj, obj in enumerate(imo["data_files"]):
        if obj["name"] == "instrument_info":
            # This is a telescope
            name = obj["metadata"]["name"]
            props = {x: obj["metadata"][x] for x in tel_keys}
            tel_data[name] = props
            for ch in props["channel_names"]:
                chan_to_tel[ch] = name
        elif obj["name"] == "channel_info":
            # This is a frequency
            name = obj["metadata"]["channel"]
            props = {x: obj["metadata"][x] for x in chan_keys}
            chan_data[name] = props
        elif obj["name"] == "scanning_parameters":
            scan_data = {x: obj["metadata"][x] for x in scan_keys}

    # Now go through detectors
    det_keys = [
        "wafer",
        "pixel",
        "pixtype",
        "channel",
        "squid",
        "psd_dac_ukcmb^2",
        "fwhm_arcmin",
        "ellipticity",
        "bandcenter_ghz",
        "bandwidth_ghz",
        "sampling_rate_hz",
        "net_ukrts",
        "pol_sensitivity_ukarcmin",
        "fknee_mhz",
        "fmin_hz",
        "alpha",
        "pol",
        "orient",
        "quat",
    ]
    obs_det_data = dict()
    wafer_to_chan = dict()
    for iobj, obj in enumerate(imo["data_files"]):
        if obj["name"] == "detector_info":
            # This is a detector
            name = obj["metadata"]["name"]
            if det_pat.match(name) is None:
                continue
            props = {x: obj["metadata"][x] for x in det_keys}
            pixstr = f"{props['pixel']:03d}"
            if pixel_pat.match(pixstr) is None:
                continue
            if wafer_pat.match(props["wafer"]) is None:
                continue
            if chan_pat.match(props["channel"]) is None:
                continue
            tel = chan_to_tel[props["channel"]]
            if tel_pat.match(tel) is None:
                continue
            if props["wafer"] not in wafer_to_chan:
                wafer_to_chan[props["wafer"]] = props["channel"]
            if wafer_obs:
                obs_key = props["wafer"]
            else:
                obs_key = props["channel"]
            if obs_key not in obs_det_data:
                obs_det_data[obs_key] = dict()
            obs_det_data[obs_key][name] = props

    # Now that we have all properties for our selected detectors, build
    # a table for each TOAST focalplane we will use.

    obs_tel = list()

    for oname, odets in obs_det_data.items():
        detlist = list(sorted(odets.keys()))
        n_det = len(detlist)

        if n_det == 0:
            msg = f"Telescope {oname} has no detectors"
            raise RuntimeError(msg)

        # Now go through the detectors polarization / orientation and set the gamma
        # and psi_pol angles...

        rate = odets[detlist[0]]["sampling_rate_hz"] * u.Hz

        if wafer_obs:
            chan = wafer_to_chan[oname]
            tel = chan_to_tel[chan]
            tname = f"{tel}_{chan}_{oname}"
        else:
            chan = oname
            tel = chan_to_tel[chan]
            tname = f"{tel}_{chan}"
        
        det_table = QTable(
            [
                Column(name="name", data=detlist),
                Column(name="wafer", data=[odets[x]["wafer"] for x in detlist]),
                Column(name="pixel", data=[odets[x]["pixel"] for x in detlist]),
                Column(name="pixtype", data=[odets[x]["pixtype"] for x in detlist]),
                Column(name="squid", data=[odets[x]["squid"] for x in detlist]),
                Column(name="quat", data=[odets[x]["quat"] for x in detlist]),
                Column(name="gamma", length=n_det, unit=u.rad),
                Column(name="psi_pol", length=n_det, unit=u.rad),
                Column(name="pol_leakage", data=np.zeros(n_det, dtype=np.float32)),
                Column(name="fwhm", length=n_det, unit=u.arcmin),
                Column(
                    name="ellipticity", data=[odets[x]["ellipticity"] for x in detlist]
                ),
                Column(name="bandcenter", length=n_det, unit=u.GHz),
                Column(name="bandwidth", length=n_det, unit=u.GHz),
                Column(name="psd_net", length=n_det, unit=(u.K * np.sqrt(1.0 * u.second))),
                Column(name="psd_fknee", length=n_det, unit=u.Hz),
                Column(name="psd_fmin", length=n_det, unit=u.Hz),
                Column(name="psd_alpha", data=[odets[x]["alpha"] for x in detlist]),
                Column(name="pol", data=[odets[x]["pol"] for x in detlist]),
                Column(name="orient", data=[odets[x]["orient"] for x in detlist]),
                Column(name="psd_dac",length=n_det, unit=u.K**2),
                Column(name="pol_sensitivity", length=n_det, unit=u.K / u.arcmin),
            ]
        )
        for idet, det in enumerate(detlist):
            det_table[idet]["gamma"] = 0.0 * u.rad
            det_table[idet]["psi_pol"] = 0.0 * u.rad
            det_table[idet]["fwhm"] = odets[det]["fwhm_arcmin"] * u.arcmin
            det_table[idet]["bandcenter"] = odets[det]["bandcenter_ghz"] * u.GHz
            det_table[idet]["bandwidth"] = odets[det]["bandwidth_ghz"] * u.GHz
            det_table[idet]["psd_net"] = odets[det]["net_ukrts"] * 1.0e-6 * (u.K * np.sqrt(1.0 * u.second))
            det_table[idet]["psd_fknee"] = odets[det]["fknee_mhz"] * 1.0e-3 * u.Hz
            det_table[idet]["psd_fmin"] = odets[det]["fmin_hz"] * u.Hz
            det_table[idet]["psd_dac"] = odets[det]["psd_dac_ukcmb^2"] * 1.0e-12 * (u.K**2)
            det_table[idet]["pol_sensitivity"] = odets[det]["pol_sensitivity_ukarcmin"] * 1.0e-6 * (u.K / u.arcmin)

        site = LitebirdSite()
        focalplane = Focalplane(
            detector_data=det_table,
            sample_rate=rate,
        )
        
        tele = Telescope(tname, focalplane=focalplane, site=site)
        tele.hwp_rpm = tel_data[tel]["hwp_rpm"]
        tele.sampling_rate = u.Quantity(tel_data[tel]["sampling_rate_hz"], unit=u.Hz)
        tele.spin_boresight_angle = u.Quantity(tel_data[tel]["spin_boresight_angle_deg"], unit=u.degree)
        tele.boresight_rotangle = u.Quantity(tel_data[tel]["boresight_rotangle_deg"], unit=u.degree)
        tele.spin_rotangle = u.Quantity(tel_data[tel]["spin_rotangle_deg"], unit=u.degree)
        obs_tel.append(tele)
    return scan_data, obs_tel


class LitebirdSchedule(SatelliteSchedule):

    def __init__(
        self, 
        scan_props, 
        mission_start, 
        observation_time, 
        num_obs=None,
    ):
        log = Logger.get()
        spin_rate_rpm = scan_props["spin_rate_rpm"]
        spin_period = (1.0 / spin_rate_rpm) * u.minute
        duty_cycle = scan_props["observation_duty_cycle"]
        duration = scan_props["mission_duration_year"] * u.year
        precession_period = scan_props["precession_period_min"] * u.minute

        obs_seconds = observation_time.to_value(u.second)
        gap_seconds = (1.0 - duty_cycle) * obs_seconds / duty_cycle
        total_seconds = obs_seconds + gap_seconds

        if num_obs is None:
            # Simulating the full mission
            num_obs = int(duration.to_value(u.second) / total_seconds)
        else:
            if num_obs * total_seconds > duration.to_value(u.second):
                new_num_obs = int(duration.to_value(u.second) / total_seconds)
                msg = f"Simulating {num_obs} observations of {total_seconds:0.2e}"
                msg += f" seconds each exceeds mission duration.  "
                msg += f"Using {new_num_obs} instead"
                num_obs = new_num_obs
                log.warning(msg)
        
        if mission_start.tzinfo is None:
            msg = f"Mission start time '{mission_start}' is not timezone-aware.  Assuming UTC."
            log.warning(msg)
            mission_start = mission_start.replace(tzinfo=timezone.utc)

        obs = timedelta(seconds=obs_seconds)
        gap = timedelta(seconds=gap_seconds)
        epsilon = timedelta(seconds=0)
        if gap_seconds == 0:
            # If there is no gap, we add a tiny break (much less than one sample for any
            # reasonable experiment) so that the start time of one observation is never
            # identical to the stop time of the previous one.
            epsilon = timedelta(microseconds=2)

        total = obs + gap

        scans = list()
        for sc in range(num_obs):
            start = sc * total + mission_start
            stop = start + obs - epsilon
            name = "{}{:06d}_{}".format("LB_", sc, start.isoformat(timespec="minutes"))
            scans.append(
                SatelliteScan(
                    name=name,
                    start=start,
                    stop=stop,
                    prec_period=precession_period,
                    spin_period=spin_period,
                )
            )
        super().__init__(
            scans=scans,
            site_name="LiteBIRD",
            telescope_name="LiteBIRD",
        )