# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:27:10 2021

@author: Libra


Deprecated as of 2021-03-12

"""

import os
import csv

import h5py
import numpy as np
from ..per_sec_stats import align

### all brain region entry point

def prepare_data_sync():
    dpath = align.get_root_path()
    syncdata = []
    for path in align.traverse(dpath):
        print(path)
        SU_ids = []
        trials = None
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
                    dset = ffr["Trials"]
                    trials = np.array(dset, dtype="double").T
                done_read = True
            except OSError:
                print("h5py read error handled")
        if trials is None:
            continue
        suid_reg = []
        with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
            l = list(csv.reader(csvfile))[1:]
            suid_reg = [list(i) for i in zip(*l)]

        (perf_desc, perf_code, welltrain_window, correct_resp) = align.judgePerformance(trials)

        if perf_code != 3:
            continue
        ### TODO: export per su per file alignment file

        for one_su in SU_ids.flatten():
            # breakpoint()
            # print(path+','+one_su.astype(np.str))
            syncdata.append([path, one_su])

    with open("su_list.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in syncdata:
            cwriter.writerow(row)

    return syncdata