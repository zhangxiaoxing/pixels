# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 09:51:12 2021

@author: Libra
"""

import re
import os

import numpy as np

def get_miceid_date_imecno(one_probe):
    datestr=re.findall(r"^\d{6,8}(?=\D)|(?<=\D)\d{6,8}(?=\D)",one_probe)[0]
    if len(datestr)==6:
        datestr="20"+datestr
    imecNo = re.findall(r"(?<=imec)\d",one_probe)
    miceId = re.findall(r"^\d{2,3}(?=\D)|(?<=\D)\d{2,3}(?=\D)",one_probe)


    if datestr and imecNo and miceId:
        return (miceId[0], datestr, imecNo[0])
    else:
        print("unresolved path meta data: ", one_probe)
        input("press Enter to continue")
        return (None, None, None)


def get_bsid_duration_who(ap_meta_path):
    bs_id = 0
    time_s = 0
    with open(ap_meta_path, "r") as file:
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


def get_root_path():
    breakpoint()
    dpath = None
    if os.path.exists("/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"):
        dpath = "/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"
    elif os.path.exists("/public1/home/sc51281/neupix/DataSum/"):
        dpath = "/public1/home/sc51281/neupix/DataSum/"
    else:
        dpath = r"K:\neupix\SPKINFO"

    return dpath


def traverse(path):
    for (basepath, dirs, files) in os.walk(path):
        if "FR_All_1000.hdf5" in files:
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

    #assuming 10-col trial array
    sample_loc = 4
    test_loc = 5
    lick_loc = 6

    if trials.shape[0] >= 40:
        correctResp = np.bitwise_xor(trials[:, sample_loc] == trials[:, test_loc], trials[:, lick_loc] == 1)
        inWindow = np.zeros((trials.shape[0],), dtype="bool")
        i = 40
        while i <= trials.shape[0]:
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
                return ("passive", 0, np.zeros_like(inWindow), correctResp)
            else:
                return ("transition", 1, np.zeros_like(inWindow), correctResp)
