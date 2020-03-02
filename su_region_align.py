# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:29:29 2019

@author: Libra

This script is responsible for generating the  su_id2reg.csv which is further used by other statistics

"""

import os
import h5py
import re
import csv
import numpy as np
import pandas as pd
import selectivity as zpy

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


def matchDepth(depth, depthL, date, mice_id, imec_no, unlabeledRecord):
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
    if re.match("CA[13]", r):
        return r
    if re.match("([A-Za-z-]+)[1-6/]{0,3}[ab]{0,1}", r):
        g = re.match("([A-Za-z-]+)[1-6/]{0,3}[ab]{0,1}", r)
        return g.group(1)
    if r.startswith("COAp"):
        return "COAp"
    else:
        return r


def getRegionList():
    site_file = r"D:\neupix\meta\NP tracks validatedFeb26.csv"
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
    return regionL


### Walk through
if __name__ == "__main__":
    regionL = getRegionList()
    unlabeledRecord = []

    for path in zpy.traverse(r"K:\neupix\DataSum"):

        # if os.path.isfile(os.path.join(path, 'su_id2reg.csv')):
        #     continue
        if not os.path.isfile(os.path.join(path, "FR_All.hdf5")):
            print("missing one FR_All file, path: ", path)
            continue

        (bs_id, time_s, who) = get_bsid_duration_who(path)
        (mice_id, date, imec_no) = zpy.get_miceid_date_imecno(path)
        depthL = getTrackRegion(regionL, mice_id, date, imec_no, who)
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
            su_region_corr.append([one_su, reg])

        tbl = pd.DataFrame(su_region_corr, columns=['index', 'region']).set_index('index')['region']
        tbl.to_csv(os.path.join(path, 'su_id2reg.csv'), header=True)

        with open("unlabeled.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(unlabeledRecord)
