# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 16:08:14 2021

@author: Libra
"""

import os
import re
import csv

import numpy as np
import h5py

from per_sec_stats.per_sec_data import per_sec_data
from per_sec_stats.get_stats import get_stats
import align.fileutil as futil

def prepare_data(delay=6,debug=False):
    sel_list = []
    wrsp_list =[]
    cluster_id_list = []
    reg_list = []
    folder_list = []
    dpath = futil.get_root_path()
    counter = 0
    for path in sorted(futil.traverse(dpath)):
        print(path)
        with h5py.File(os.path.join(path, "FR_All_1000.hdf5"), "r") as ffr:
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

        if (trials is None) or np.sum(trials[:,8])<80:
            continue

        if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
            continue
        with open(os.path.join(path, "su_id2reg.csv"),'r') as csvfile:
            reg_l = list(csv.reader(csvfile))
            reg_l=reg_l[1:] # discard header

        dict_stats=get_stats(trial_FR, trials, delay=delay, debug=debug)

        if debug:
            SU_ids=SU_ids[:,:20]
            reg_l=reg_l[:20]

        sel_list.extend(dict_stats['per_sec_selectivity']) #(n(Delay), n(SU))
        wrsp_list.extend(dict_stats['per_sec_wrs_p']) #(n(Delay), n(SU))
        cluster_id_list.extend(SU_ids[0]) # (n(SU),)
        folder_list.extend([re.search(r"(?<=SPKINFO\\)(.*)", path)[1]] * SU_ids.shape[1]) # (n(SU),)
        reg_list.extend(reg_l) # (n(SU),)
        counter += 1
        if debug and counter>2:
            break

    return {'cid':cluster_id_list,
            'sel':sel_list,
            'wrsp':wrsp_list,
            'reg':reg_list,
            'folder':folder_list,
            'delay':delay}