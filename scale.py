# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 13:48:45 2020

@author: Libra
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import selectivity as zpy
import h5py
from per_region_roc import Auc_stats
from GLM_delay_stats import GLM_stats


def process_one_path(path):
    # %% additional metrics
    print(path)
    unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), sep="\t")
    SU_ids = []
    trial_FR = None
    trials = None
    if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
        return
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
        return

    (perf_desc, perf_code, welltrain_window, correct_resp) = zpy.judgePerformance(trials)

    pct_stats = GLM_stats()
    auc_stats = Auc_stats()

    pct_stats.processGLMStats(trial_FR, trials, welltrain_window, correct_resp)
    auc_stats.process_auc_stats(trial_FR, trials, welltrain_window, correct_resp)

    pct_features = pct_stats.getFeatures()
    auc_features = auc_stats.get_features()

    ### TODO: features to coordinates

    # %% original scale ref

    su_ids = np.load(os.path.join(path, 'spike_clusters.npy'))
    spkTS = np.load(os.path.join(path, 'spike_times.npy'))

    su_ids = su_ids[np.squeeze((spkTS >= 30000 * 3600) & (spkTS < 30000 * (3600 + 2400)))]
    spkTS = spkTS[np.squeeze((spkTS >= 30000 * 3600) & (spkTS < 30000 * (3600 + 2400)))]

    id_set = np.unique(unitInfo['id'])
    fh = plt.figure(figsize=(4.8, 9.6))
    depth_list = []
    for one_id in id_set:
        depthDF = unitInfo.loc[unitInfo['id'] == one_id, 'depth']
        if depthDF.empty:
            print(one_id)
            continue
        depth = depthDF.iat[0]
        depth_list.append([one_id, depth])
        xs = spkTS[su_ids == one_id] / 30000
        ys = np.ones_like(xs) * depth
        plt.scatter(xs, ys, s=8, c='k', marker='|', edgecolors='none', alpha=0.002)
    fh.savefig(os.path.join(path, 'scale.png'), dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    DEBUGGING = True

    if DEBUGGING:
        path = r'K:\neupix\DataSum\191015-DPA-Learning2_29_g0\191015-DPA-Learning2_29_g0_imec0_cleaned'
        process_one_path(path)

    else:
        from multiprocessing import Pool

        curr_pool = Pool(processes=12)

        dpath = None
        if os.path.exists("/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"):
            dpath = "/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"
        elif os.path.exists(r"K:\neupix\DataSum"):
            dpath = r"K:\neupix\DataSum"
        else:
            dpath = r"D:\neupix\DataSum"

        all_proc = []
        for path in zpy.traverse(dpath):
            print(path)
            all_proc.append(curr_pool.apply_async(process_one_path, args=(path,)))

        for one_proc in all_proc:
            one_proc.get()
            print('get')
