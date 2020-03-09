"""
Created on Wed Mar 05 00:14:22 2020

@author: Libra

Generate data set for CQ's time permute transient coding algorithm (transient\transient.m)
Only produce intermediate data

"""

import numpy as np
import h5py
import sys
import os

def create_6s_dataset():
    fstr = np.load("ctd.npz", allow_pickle=True)
    features_per_su = fstr["features_per_su"].tolist()
    reg_list = fstr["reg_list"].tolist()

    with h5py.File("su_trials_fr_6.hdf5", "w") as fw:
        fw.create_dataset('count',data=len(features_per_su))
        for suid in np.arange(0, len(features_per_su)):
            S1_raw = features_per_su[suid]['S1_6']
            S2_raw = features_per_su[suid]['S2_6']

            S1 = np.array([np.mean(S1_raw[bin_idx:bin_idx + 4, :], axis=0) for bin_idx in np.arange(16, 40, 4)])
            S2 = np.array([np.mean(S2_raw[bin_idx:bin_idx + 4, :], axis=0) for bin_idx in np.arange(16, 40, 4)])

            grp = fw.create_group(suid.__str__())
            grp.create_dataset('S1', data=S1)
            grp.create_dataset('S2', data=S2)
#

def create_3s_dataset():
    fstr = np.load("ctd.npz", allow_pickle=True)
    features_per_su = fstr["features_per_su"].tolist()
    reg_list = fstr["reg_list"].tolist()

    with h5py.File("su_trials_fr_3.hdf5", "w") as fw:
        fw.create_dataset('count',data=len(features_per_su))
        for suid in np.arange(0, len(features_per_su)):
            S1_raw = features_per_su[suid]['S1_3']
            S2_raw = features_per_su[suid]['S2_3']

            S1 = np.array([np.mean(S1_raw[bin_idx:bin_idx + 4, :], axis=0) for bin_idx in np.arange(16, 28, 4)])
            S2 = np.array([np.mean(S2_raw[bin_idx:bin_idx + 4, :], axis=0) for bin_idx in np.arange(16, 28, 4)])

            grp = fw.create_group(suid.__str__())
            grp.create_dataset('S1', data=S1)
            grp.create_dataset('S2', data=S2)


if __name__=="__main__":

    if not os.path.isfile("per_sec_sel.npz"):
        print('missing data file!')
        sys.exit(0)
    fstr = np.load('per_sec_sel.npz')
    per_sec_sel_arr = fstr['per_sec_sel_arr']
    non_sel_mod_arr = fstr['non_sel_mod_arr']
    # perfS1_arr = fstr['perfS1_arr']
    # perfS2_arr = fstr['perfS2_arr']
    # reg_arr = fstr['reg_arr']
    transient6=None
    transient3=None
    with h5py.File(os.path.join('transient','CQ_transient.hdf5'), 'r') as fr:
        transient6 = np.array(fr["transient6"])
        transient3 = np.array(fr["transient3"])