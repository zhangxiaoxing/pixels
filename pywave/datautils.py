# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:03:16 2021

@author: Libra
"""

import os
import platform
import h5py
import numpy as np

if __package__:
    import align.fileutil as fileutil
else:
    import sys
    sys.path.append(os.path.dirname(__file__)+'/../align/fileutil')
    import fileutil as fileutil




def traverse_data():
    if platform.system() == "Linux":
        homedir = r'/home/zx/neupix/SPKINFO'
    elif platform.system() == "Windows":
        homedir = r'k:\neupix\SPKINFO'
    for (basepath, dirs, files) in os.walk(homedir):
        if "FR_All_ 250.hdf5" in files:
            yield basepath


def get_dataset(correct_error='correct'):
    fr_per_su = []
    # avails = []
    # reg_list = []

    for path in traverse_data():
        print(path)
        # if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
        #     continue

        with h5py.File(os.path.join(path, "FR_All_ 250.hdf5"), "r") as ffr:
            SU_ids = np.array(ffr["SU_id"],dtype=np.uint16)
            trial_FR = np.array(ffr["FR_All"], dtype=np.int32)
            trials = np.array(ffr["Trials"], dtype=np.int64)

        # if trials is None:
        #     continue
        (_perf_desc, perf_code, welltrain_window, correct_resp,
         ) = fileutil.judgePerformance(trials.transpose())
        # if np.sum(trials[8,:])<40:
        if perf_code != 3:
            continue

        (s1_3_sel, s1_6_sel, s2_3_sel, s2_6_sel) = get_trials(correct_error, trials)

        onesu = {'S1_3': trial_FR[:, :, s1_3_sel],
                 'S1_6': trial_FR[:, :, s1_6_sel],
                 'S2_3': trial_FR[:, :, s2_3_sel],
                 'S2_6': trial_FR[:, :, s2_6_sel],
                 'suid':SU_ids,}

        fr_per_su.append(onesu)

    return (fr_per_su)


def get_meta(delay=6):
    homedir=getHomedir()
    with h5py.File(os.path.join(homedir, f'transient_{delay}.hdf5'), "r") as ffr:
        allcid = np.array(ffr["cluster_id"],dtype=np.uint16)
        mem_type = np.array(ffr["mem_type"], dtype=np.int32)

    meta = {'allcid': allcid,
             'mem_type': mem_type,
             }
    return meta

def getHomedir(type='sums'):
    if sys.platform=='win32':
        if type=='sums':
            homedir = os.path.join('K:','code','per_sec');
        elif type=='raw':
            homedir = os.path.join('K:','neupix','SPKINFO');

    elif sys.platform=='linux':
        if type=='sums':
            homedir = os.path.join('~','pixels','per_sec');
        elif type=='raw':
            homedir = os.path.join('~','neupix','SPKINFO');

    return homedir

    # suid_reg = []
    # TODO brain region
    # with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
    #     l = list(csv.reader(csvfile))[1:]
    #     suid_reg = [list(i) for i in zip(*l)]


def get_trials(correct_error, trials):
    if correct_error == 'correct':
        s1_3_sel = np.all([trials[4, :] == 4, trials[7, :] ==
                           3, trials[8, :], trials[9, :]], axis=0)
        s1_6_sel = np.all([trials[4, :] == 4, trials[7, :] ==
                           6, trials[8, :], trials[9, :]], axis=0)
        s2_3_sel = np.all([trials[4, :] == 8, trials[7, :] ==
                           3, trials[8, :], trials[9, :]], axis=0)
        s2_6_sel = np.all([trials[4, :] == 8, trials[7, :] ==
                           6, trials[8, :], trials[9, :]], axis=0)
    else:
        s1_3_sel = np.all([trials[4, :] == 4, trials[7, :] ==
                           3, trials[9, :] == 0], axis=0)
        s1_6_sel = np.all([trials[4, :] == 4, trials[7, :] ==
                           6,  trials[9, :] == 0], axis=0)
        s2_3_sel = np.all([trials[4, :] == 8, trials[7, :] ==
                           3,  trials[9, :] == 0], axis=0)
        s2_6_sel = np.all([trials[4, :] == 8, trials[7, :] ==
                           6,  trials[9, :] == 0], axis=0)
    return (s1_3_sel, s1_6_sel, s2_3_sel, s2_6_sel)