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

from per_sec_stats.get_stats import get_stats
import align.fileutil as futil

def prepare_data(delay=6,debug=False):
    '''
    Assign per SU selectivity and brain-region-tree for further analysis


    Parameters
    ----------
    delay : TYPE, optional
        WM task delay duration. The default is 6.
    debug : TYPE, optional
        Developer only. The default is False.

    Returns
    -------
    dict
        Per SU selectivity and brain-region-tree.
    error_files : TYPE
        Sub-criteria performance and incomplete data. Developer only.

    '''

    error_files=[]
    sel_list = []
    wrsp_list =[]
    trialconts_list=[]
    cluster_id_list = []
    wf_good_list = []
    reg_list = []
    folder_list = []
    dpath = futil.get_root_path()
    counter = 0
    for path in sorted(futil.traverse(dpath)):
        print(path)
        with h5py.File(os.path.join(path, "FR_All_1000.hdf5"), "r") as ffr: # Trial-aligned firing rate
            # print(list(ffr.keys()))
            if not "SU_id" in ffr.keys():
                done_read = True
                # print("missing su_id key in path ", path)
                error_files.append(['Missing SU_id',path])
                continue
            SU_ids = np.array(ffr["SU_id"], dtype="uint16")
            trial_FR = np.array(ffr["FR_All"], dtype="double")
            trials = np.array(ffr["Trials"], dtype="double").T
            WF_good = np.array(ffr["WF_good"], dtype="int8") # waveform criteria

        if np.sum(trials[:,8])<40: #Well-trained criteria
            error_files.append(['Lack trials',path])
            continue

        if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
            error_files.append(['Missing localization',path]) #Localization failure, prob. missing IHC data
            continue

        with open(os.path.join(path, "su_id2reg.csv"),'r') as csvfile:
            reg_l = list(csv.reader(csvfile))
            reg_l=reg_l[1:] # discard header

        dict_stats=get_stats(trial_FR, trials, delay=delay, debug=debug) #FR statistics

        if debug:
            SU_ids=SU_ids[:,:20]
            reg_l=reg_l[:20]
            WF_good=WF_good[:,:20]

        reg_suids=[np.uint16(np.float64(x[0])) for x in reg_l]
        if not np.all([reg_suids[x]==SU_ids[0][x] for x in range(len(reg_suids))]): # Assuming identical sort to save time
            breakpoint()
        # (map->) reduce to single dataset
        sel_list.extend(dict_stats['per_sec_selectivity']) #(n(SU), n(Delay_bin))
        wrsp_list.extend(dict_stats['per_sec_wrs_p']) #(n(SU), n(Delay_bin))
        trialconts_list.extend(dict_stats['trial_counts']) #(n(SU), [s1, s2])
        cluster_id_list.extend(SU_ids[0]) # (n(SU),)
        wf_good_list.extend(WF_good[0]) # (n(SU),)
        folder_list.extend([re.search(r"(?<=SPKINFO\\)(.*)", path)[1]] * SU_ids.shape[1]) # (n(SU),)
        reg_list.extend(reg_l) # (n(SU),)
        counter += 1
        # if debug and counter>2:
        #     break

    return ({'cid':cluster_id_list,
             'WF_good':wf_good_list,
            'sel':sel_list,
            'wrsp':wrsp_list,
            'reg':reg_list,
            'folder':folder_list,
            'trial_counts':trialconts_list,
            'delay':delay},
            error_files)