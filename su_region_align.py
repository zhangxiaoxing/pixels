# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:29:29 2019

@author: Libra

This script is responsible for generating the  su_id2reg.csv which is further used by other statistics
first column is cluster id, second is current alignment, third is old alignment


"""

import os
import h5py
import sys
import re
import csv
import numpy as np
import pandas as pd
import selectivity as zpy


def traverse(path):
    for (basepath, dirs, files) in os.walk(path):
        if "cluster_info.tsv" in files:
            yield basepath


def judgePerformance(trials, criteria=75):
    """

    Parameters
    ----------
    trials : TYPE
        behaviral trials array.

    Returns
    -------
    str
        readable description of the learning stage.
    int
        single integer code for the learning stage.
    trials
        if well-trained, the trials in the engaging window, all trials otherwise.

    """
    sample_loc = 4 if trials.shape[1] > 6 else 2
    test_loc = 5 if trials.shape[1] > 6 else 3
    lick_loc = 6 if trials.shape[1] > 6 else 4

    if trials.shape[0] >= 40:
        correctResp = np.bitwise_xor(trials[:, sample_loc] == trials[:, test_loc], trials[:, lick_loc] == 1)
        inWindow = np.zeros((trials.shape[0],), dtype="bool")
        i = 40
        while i < trials.shape[0]:
            #            if np.sum(correctResp[i-40:i])>=32:
            if np.sum(correctResp[i - 40: i]) >= criteria * 40 / 100:
                inWindow[i - 40: i] = 1
            i += 1
        if np.sum(inWindow) >= 40:  # Well Trained
            return ("wellTrained", 3, inWindow, correctResp)

        else:
            inWindow = np.zeros((trials.shape[0],), dtype="bool")
            licks = trials[:, lick_loc] == 1
            i = 40
            while i < trials.shape[0]:
                if np.sum(licks[i - 40: i]) >= 16:  # Learning
                    inWindow[i - 40: i] = 1
                i += 1
            if np.sum(inWindow) >= 40:
                return ("learning", 2, inWindow, correctResp)
            elif np.sum(licks) <= trials.shape[0] // 10:  # Passive
                return ("passive", 0, np.ones_like(inWindow), correctResp)
            else:
                return ("transition", 1, np.ones_like(inWindow), correctResp)


def get_root_path():
    dpath = None
    if os.path.exists("/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"):
        dpath = "/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"
    elif os.path.exists(r"K:\neupix\DataSum"):
        dpath = r"K:\neupix\DataSum"
    elif os.path.exists("/public1/home/sc51281/neupix/DataSum/"):
        dpath = "/public1/home/sc51281/neupix/DataSum/"
    else:
        dpath = r"D:\neupix\DataSum"

    return dpath


def get_bsid_duration_who(path):
    bs_id = 0
    time_s = 0
    files = os.listdir(path)
    for f in files:
        if f.endswith("ap.meta"):
            with open(os.path.join(path, f), "r") as file:
                for line in file:
                    bs_id_grps = re.match("imDatBsc_sn=(\\d{1,3})", line)
                    if bs_id_grps:
                        bs_id = bs_id_grps.group(1)
                    fs_grps = re.match("fileSizeBytes=(\\d+)", line)
                    if fs_grps:
                        time_s = int(fs_grps.group(1)) / 385 / 2 / 30000

    if bs_id == "350":
        who_did = "ZHA"
    elif bs_id == "142":
        who_did = "HEM"
    else:
        who_did = "UNKNOWN"
        print("Unknown BS id!")
        input("press Enter to continue")

    return (bs_id, time_s, who_did)


def imecNo2side(who_did, date, imecNo, mid):
    if date == "191130" and mid == "49":
        return "L"

    if date == "191101" and mid == "26":
        return "R"

    if who_did == "HEM" and (int(date)) >= 191028:
        if imecNo == "1":
            return "R"
        elif imecNo == "0":
            return "L"
    else:
        if imecNo == "1":
            return "L"
        elif imecNo == "0":
            return "R"

    print("Error parsing imec No")
    return "X"


def matchDepth(depth, depthL, date, mice_id, imec_no, unlabeledRecord=None):
    if not depthL.empty:
        label = depthL.loc[
            (depthL["distance2tipLow"] <= depth) & (depth < depthL["distance2tipHigh"]),
            ["acronym"],
        ]
        if len(label.index) == 1:
            return label.iloc[0, 0]
    if unlabeledRecord is not None:
        unlabeledRecord.append([date, mice_id, imec_no, depth])
    return "Unlabeled"


def getTrackRegion(regionL, mice_id, date, imecNo, who_did):
    depthL = regionL.loc[
        (regionL["mouse_id"] == mice_id)
        & (regionL["implanting_date"] == "20" + date)
        & (regionL["side"] == imecNo2side(who_did, date, imecNo, mice_id)),
        ["acronym", "distance2tipLow", "distance2tipHigh"],
    ]
    return depthL


def combineSubRegion(r):
    if re.match("CA[123]", r):
        return r
    if r.startswith("COAp"):
        return "COAp"
    if r.startswith("DG-"):
        return "DG"
    if r.startswith("LGd-"):
        return "LGd"
    if r.startswith("SSp-"):
        return "SSp"
    if re.match("([A-Za-z-]+)[1-6/]{0,3}[ab]{0,1}", r):
        g = re.match("([A-Za-z-]+)[1-6/]{0,3}[ab]{0,1}", r)
        return g.group(1)
    else:
        return r


def getRegionList(cmp=False):
    site_file = r"K:\neupix\meta\NP tracks revised 0319.csv"
    regionL = pd.read_csv(site_file).astype(
        {"mouse_id": "str", "implanting_date": "str"}
    )[
        [
            "mouse_id",
            "implanting_date",
            "side",
            "acronym",
            "distance2tipHigh",
            "distance2tipLow",
        ]
    ]

    regionL["acronym"] = regionL["acronym"].apply(combineSubRegion)
    if not cmp:
        return (regionL,None)
    else:
        cmp_file = r"K:\neupix\meta\NP tracks validated3.3.csv"
        cmpL = pd.read_csv(cmp_file).astype(
            {"mouse_id": "str", "implanting_date": "str"}
        )[
            [
                "mouse_id",
                "implanting_date",
                "side",
                "acronym",
                "distance2tipHigh",
                "distance2tipLow",
            ]
        ]

        cmpL["acronym"] = cmpL["acronym"].apply(combineSubRegion)
        return (regionL, cmpL)


### Walk through
if __name__ == "__main__":
    cmp=False
    identical_counter = []
    (regionL, cmpL) = getRegionList(cmp=cmp)
    unlabeledRecord = []

    for path in traverse(r"K:\neupix\DataSum"):

        # if os.path.isfile(os.path.join(path, 'su_id2reg.csv')):
        #     continue
        if not os.path.isfile(os.path.join(path, "FR_All.hdf5")):
            print("missing one FR_All file, path: ", path)
            continue

        print(path)
        (bs_id, time_s, who) = get_bsid_duration_who(path)
        (mice_id, date, imec_no) = zpy.get_miceid_date_imecno(path)
        depthL = getTrackRegion(regionL, mice_id, date, imec_no, who)
        if cmp:
            depth_cmp = getTrackRegion(cmpL, mice_id, date, imec_no, who)
        su_ids = None
        su_region_corr = []
        with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
            if not "SU_id" in ffr.keys():
                print("missing su_id key in path ", path)
                continue
            dset = ffr["SU_id"]
            su_ids = np.array(dset, dtype="uint16")[0]

        unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), sep="\t")

        ### su ids match reg in sequence
        for one_su in su_ids:
            depth = unitInfo.loc[unitInfo['id'] == one_su, ['depth']].iat[0, 0]
            reg = matchDepth(depth, depthL, date, mice_id, imec_no, unlabeledRecord)
            if cmp:
                reg_cmp = matchDepth(depth, depth_cmp, date, mice_id, imec_no, unlabeledRecord)
                su_region_corr.append([one_su, reg, reg_cmp])
                if reg == reg_cmp:
                    identical_counter.append([1, reg, reg_cmp])
                else:
                    identical_counter.append([0, reg, reg_cmp])
            else:
                su_region_corr.append([one_su, reg])
        if cmp:
            tbl = pd.DataFrame(su_region_corr, columns=['index', 'region', 'old']).set_index('index')[:]
        else:
            tbl = pd.DataFrame(su_region_corr, columns=['index', 'region']).set_index('index')[:]
        tbl.to_csv(os.path.join(path, 'su_id2reg.csv'), header=True)

    with open("unlabeled.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(unlabeledRecord)

        # np.sum(list(zip(*identical_counter))[0])
    if cmp:
        print(np.sum(list(zip(*identical_counter))[0]) / len(identical_counter))
