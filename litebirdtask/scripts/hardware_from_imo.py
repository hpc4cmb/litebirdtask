# Copyright (c) 2015-2021 LiteBIRD Collaboration.
# Full license can be found in the top level "LICENSE" file.
"""
Plot a hardware model.
"""

import argparse
import json
import numpy as np
import toml

def get_wafer_dic (tele):

    wnp = {
            "LF1": 9,
            "LF2": 36,
            "MF1": 61,
            "MF2": 61,
            "HF1": 127,
            "HF2": 127,
            "HF3": 169,
        }

    pix2w = {
            "LF1": ["LP1", "LP2"],
            "LF2": ["LP3", "LP4"],
            "MF1": ["MP1"],
            "MF2": ["MP2"],
            "HF1": ["HP1"],
            "HF2": ["HP2"],
            "HF3": ["HP3"],
        }

    toast_keys = ['pixels', 'bands', 'telescopes', 'wafers', 'detectors']
    pixmm = {
            "LP1": 32.0,
            "LP2": 32.0,
            "LP3": 16.0,
            "LP4": 16.0,
            "MP1": 11.6,
            "MP2": 11.6,
            "HP1": 6.6,
            "HP2": 6.6,
            "HP3": 5.7,
        }
    pbd = {
            "LP1": ["L1-040", "L1-060", "L1-078"],
            "LP2": ["L2-050", "L2-068", "L2-089"],
            "LP3": ["L3-068", "L3-089", "L3-119"],
            "LP4": ["L4-078", "L4-100", "L4-140"],
            "MP1": ["M1-100", "M1-140", "M1-195"],
            "MP2": ["M2-119", "M2-166"],
            "HP1": ["H1-195", "H1-280"],
            "HP2": ["H2-235", "H2-337"],
            "HP3": ["H3-402"],
        }


    windx = 0
    wafers={}
    wf =  {}

    if tele[0]== "L":
        for wt, wnum, hd in zip(
                ["LF1", "LF2", "LF2", "LF1", "LF1", "LF2", "LF2", "LF1"],
                [0, 1, 2, 3, 4, 5, 6, 7],
                [0, 0, 1, 1, 1, 1, 0, 0],
            ):
            wn = "L{:02d}".format(wnum)
            wf["type"] = wt
            wf["npixel"] = wnp[wt]
            wf["pixels"] = pix2w[wt]
            wf["handed"] = hd
            wf["position"] = windx
            wf["telescope"] = "LFT"
            wafers[wn] = wf
            windx += 1

    elif tele[0]== "M":
        for wt, wnum in zip(
                ["MF1", "MF2", "MF2", "MF1", "MF2", "MF2", "MF1"], [3, 1, 0, 2, 5, 6, 4]
            ):
                wn = "M{:02d}".format(wnum)
                wf["type"] = wt
                wf["npixel"] = wnp[wt]
                wf["pixels"] = pix2w[wt]
                wf["position"] = windx
                wf["telescope"] = "MFT"
                wafers[wn] = wf
                windx += 1
    elif tele[0] =="H":
        for wt, wnum in zip(["HF1", "HF2", "HF3"], [0, 1, 2]):
            wn = "H{:02d}".format(wnum)
            wf["type"] = wt
            # wf["rhombusgap"] = wrhombgap[wt]
            wf["npixel"] = wnp[wt]
            wf["pixels"] = pix2w[wt]
            wf["position"] = windx
            wf["telescope"] = "HFT"
            wafers[wn] = wf
            windx += 1
    return wafers



def main():
    parser = argparse.ArgumentParser(
        description="This program reads the IMO `json` file downloaded from the \
        Instrument MOdel repo and makes it compatible with toast3 focalplane.\
        It produces hardware file per frequency channel. ",

        usage="lbt_hardware_from_imo  ",
    )
    parser.add_argument(
    "--hardware", required=True, default=None, help="Input json file"
    )
    args = parser.parse_args()

    print("Loading hardware file {}...".format(file), flush=True)
    with open(args.hardware) as json_file:
        data = json.load(json_file)


    # Loop over the entries of json file, populating the new dictionary per freq. channel

    for i in range(len(data ['data_files']) ) :

        if data['data_files'][i]['name'] == "instrument_info":
            instrument_dic = data['data_files'][i]['metadata']
            toast_wf =get_wafer_dic( instrument_dic['name'])
            continue
        elif data['data_files'][i]['name'] == "channel_info":

            channel_dic =data['data_files'][i] ['metadata']
            newhw={k:{}  for k in toast_keys }
            iold =i
            continue
        elif data['data_files'][i]['name'] == "detector_info" :
            det_dic=(data['data_files'][i] ['metadata']  )
        else:
            break

        newhw ['telescopes']  [ instrument_dic['name'] ] ={'wafers' : instrument_dic['wafers'],
                                                          'platescale':instrument_dic['platescale_deg_mm'],
                                                          'waferspace' : instrument_dic['waferspace_mm']
                                                              }
        newhw ['pixels']   [det_dic['pixtype']]=  {
                                                 'bands' :pbd[det_dic['pixtype']] ,
                                                'sizemm':pixmm[det_dic['pixtype']]
                                                }


        newhw ['wafers'][det_dic['wafer']] =   {
                                            'type':toast_wf[det_dic['wafer']]['type'],
                                            'npixel':np.int_(channel_dic['number_of_detectors']/2),
                                            'pixels':[det_dic['pixtype']],
                                            'position':toast_wf[det_dic['wafer']]['position'],
                                            'telescope': instrument_dic['name']
                                             }

        newhw ['bands'] [ channel_dic['channel']]= {
            'center': channel_dic['bandcenter_ghz'],
            'bandwidth': channel_dic['bandwidth_ghz'],
            'bandpass':'',
            'NET': channel_dic['net_detector_ukrts'],
            'fwhm': channel_dic['fwhm_arcmin'],
            'fknee': channel_dic['fknee_mhz'],
            'fmin': channel_dic['fmin_hz'],
            'alpha': channel_dic['alpha']}


        newhw['detectors'][det_dic['name']] = { 'wafer': det_dic['wafer'],
                              'pixel': f"{det_dic['pixel']:03d}",
                              'pixtype':det_dic[ "pixtype" ],
                              'band': det_dic[ "channel"],
                              'fwhm': det_dic[ "fwhm_arcmin"],
                              'pol':det_dic[ "pol" ],
                              'orient':det_dic[ "orient" ],
                              'quat':det_dic[ "quat" ],
                              'UID':data['data_files'][i]['uuid']
                             }
        if i - iold == channel_dic['number_of_detectors']:
            print(f"Dumping {i-iold} detectors into {channel_dic['channel']}.toml.gz"  )
            with open( f"{channel_dic['channel']}.toml.gz", "w") as f:
                toml.dump( newhw, f)

            break 

    return
