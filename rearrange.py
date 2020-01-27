# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 15:48:09 2020

@author: Libra
"""


import os
import h5py
import re
import numpy as np
import selectivity as zpy
import pandas as pd
import matplotlib.pyplot as plt

def rearrange_row(trials, trial_FR):
    row_sel_6 = np.concatenate((np.arange(16), np.arange(16, 40, 2), np.arange(40, 68)))
    row_sel_3 = np.arange(56)
    row_arr_out=[]
    for i in range(trials.shape[0]):
        if trials[i,5]<=3:
            row_arr_out.append(trial_FR[row_sel_3,i,:])
        elif trials[i,5]>=6:
            row_arr_out.append(trial_FR[row_sel_6,i,:])
        else:
            print('Error delay duration')
            breakpoint()
            
    return np.moveaxis(np.array(row_arr_out),2,0)
            
                
    
    

def rearrange_block(trials,upper_bound=0):
    if upper_bound==0:
        upper_bound=trials.shape[0]
    block_begin = 0
    all_out = np.zeros(upper_bound, dtype=np.int32)
    while block_begin < upper_bound:
        unregistered = np.full(8, True)
        for i in range(block_begin, block_begin + 8):
            if trials[i, 2] == 4 and trials[i, 3] == 4:
                if unregistered[0]:
                    all_out[i] = block_begin
                    unregistered[0] = False
                else:
                    all_out[i] = block_begin + 1
                    unregistered[4] = False
            elif trials[i, 2] == 4 and trials[i, 3] == 8:
                if unregistered[1]:
                    all_out[i] = block_begin + 2
                    unregistered[1] = False
                else:
                    all_out[i] = block_begin + 3
                    unregistered[5] = False
            elif trials[i, 2] == 8 and trials[i, 3] == 4:
                if unregistered[2]:
                    all_out[i] = block_begin + 4
                    unregistered[2] = False
                else:
                    all_out[i] = block_begin + 5
                    unregistered[6] = False
            else:  # 8, 8
                if unregistered[3]:
                    all_out[i] = block_begin + 6
                    unregistered[3] = False
                else:
                    all_out[i] = block_begin + 7
                    unregistered[7] = False
        if np.any(unregistered):
            return (False, all_out)

        block_begin += 8

    return (True, all_out)


def normalize(all_sess_arr): # SU,trials, bins
    arrmin=np.amin(all_sess_arr,axis=(1,2))
    arrspan=np.amax(all_sess_arr,axis=(1,2))-arrmin
    all_sess_arr=np.moveaxis(all_sess_arr,0,2)
    all_sess_arr-=arrmin
    all_sess_arr/=arrspan
    return np.moveaxis(all_sess_arr,2,0)
    
    
    
ntrialsCount = []
all_sess_list = []

for path in zpy.traverse("K:/neupix/DataSum/"):
    print(path)
    SU_ids = []
    trial_FR = []
    trials = []
    with h5py.File(os.path.join(path, "FR_All.hdf5")) as ffr:
        # print(list(ffr.keys()))
        if not "SU_id" in ffr.keys():
            print("missing su_id key in path ", path)
            continue
        dset = ffr["SU_id"]
        SU_ids = np.array(dset, dtype="uint16")
        dset = ffr["FR_All"]
        trial_FR = np.array(dset, dtype="double")
        dset = ffr["Trials"]
        trials = np.array(dset, dtype="double").T

    (perf_desc, perf_code, inWindow, correct_resp) = zpy.judgePerformance(trials)
    #  toReturn.extend(["wellTrained", 3, correctResp,welltrain_window])

    if perf_code != 3:
        continue
    ntrialsCount.append(trials.shape[0])
    
    if trials.shape[0] < 200:
        continue

    
    
    (reg_all, matched_index) = rearrange_block(trials,200)

    if reg_all:
        merged=rearrange_row(trials,trial_FR)
        # onesession=merged[:,matched_index,:].reshape((merged.shape[0],-1))
        onesession=merged[:,matched_index,:]
        all_sess_list.append(onesession)
        
all_sess_arr=np.concatenate(tuple(all_sess_list),axis=0)        
del all_sess_list
import gc
gc.collect()

all_sess_arr=normalize(all_sess_arr)



