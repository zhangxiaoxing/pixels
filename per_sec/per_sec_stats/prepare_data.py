# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 16:08:14 2021

@author: Libra
"""

import os
import re
import csv
import platform
import numpy as np
import h5py
import sys
from per_sec_stats.per_sec_data import per_sec_data

from align import su_region_align as align


def prepare_data(delay=6):
    curr_stats = per_sec_data()
    per_sec_sel_list = []
    non_sel_mod_list = []
    bs_sel_list = []
    perfS1_list = []
    perfS2_list = []
    raw_sel_list = []
    auc_list = []
    fr_list=[]
    wrs_p_list=[]
    # non_mod_list = []
    folder_list = []
    cluster_id_list = []
    reg_list = []
    dpath = align.get_root_path()
    for path in sorted(align.traverse(dpath)):
        print(path)
        if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
            continue
        done_read = False
        while not done_read:
            try:
                with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
                    # print(list(ffr.keys()))
                    if not "SU_id" in ffr.keys():
                        done_read = True
                        print("missing su_id key in path ", path)
                        continue
                    dset = ffr["SU_id"]
                    SU_ids = np.array(dset, dtype="uint16")
                    dset = ffr["FR_All"]
                    trial_FR = np.array(dset, dtype="double")
                    dset = ffr["Trials"]
                    trials = np.array(dset, dtype="double").T
                done_read = True
            except OSError:
                print("h5py read error handled")
        if trials is None:
            continue
        with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
            l = list(csv.reader(csvfile))[1:]
            suid_reg = [list(i) for i in zip(*l)]

        (perf_desc, perf_code, welltrain_window, correct_resp) = align.judgePerformance(trials)

        if perf_code != 3:
            continue

        curr_stats.process_select_stats(trial_FR, trials, welltrain_window, correct_resp, delay=delay)

        (per_sec_sel, non_sel_mod, perfS1, perfS2, bs_sel, raw_sel,auc,wrsp,fr) = curr_stats.getFeatures()
        per_sec_sel_list.append(per_sec_sel)
        non_sel_mod_list.append(non_sel_mod)
        bs_sel_list.append(bs_sel)
        # non_mod_list.append(non_sel_mod)
        perfS1_list.append(perfS1)
        perfS2_list.append(perfS2)
        raw_sel_list.append(raw_sel)
        auc_list.append(auc)
        fr_list.append(fr)
        wrs_p_list.append(wrsp)
        cluster_id_list.append(SU_ids)
        folder_list.append([re.search(r"(?<=DataSum\\)(.*)", path)[1]] * SU_ids.shape[1])

        for i in range(SU_ids.shape[1]):
            if suid_reg[0][i] != str(SU_ids[0, i]):
                print(f"unmatch cluster id and id_reg, {path}, {SU_ids[0, i]}")
                input("press Enter to continue")
        reg_list.append(suid_reg[1])
    return (
        per_sec_sel_list, non_sel_mod_list, perfS1_list, perfS2_list, reg_list, bs_sel_list, raw_sel_list,
        cluster_id_list, folder_list,auc_list,wrs_p_list,fr_list)