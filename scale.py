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
from GLM_delay_stats import GLM_delay_stats
from zxStats import zxStats
import su_region_align as align


# we got min=20 max=3840
def max_min_depth(path):
    unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), sep="\t")
    sessMax = max(unitInfo['depth'].values)
    sessMin = min(unitInfo['depth'].values)

    return (sessMax, sessMin)


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

    (perf_desc, perf_code, welltrain_window, correct_resp) = align.judgePerformance(trials)
    # trial_FR [68,241,157] trials[241,6] w_w,c_r [241,]

    pct_stats = GLM_delay_stats()
    auc_stats = Auc_stats()
    mod_stats = zxStats()

    pct_stats.processGLMStats(trial_FR, trials, welltrain_window, correct_resp)
    auc_stats.process_auc_stats(trial_FR, trials, welltrain_window, correct_resp)
    mod_features = mod_stats.get_normalized_fr(trial_FR, trials, welltrain_window, correct_resp)

    pct_features = pct_stats.getFeatures()  # 13 features
    auc_features = auc_stats.get_features()  # 6 features

    depth_list = [unitInfo.loc[unitInfo['id'] == id, 'depth'].values[0] for id in SU_ids[0]]

    feat_mat = np.full((192, 19), np.nan)
    mod_mat = np.full((192, 14), np.nan)
    for i in np.arange(0, 192):
        curr_depth = i * 20 + 20
        idx = [np.abs(d - curr_depth) <= 25 for d in depth_list]
        if np.count_nonzero(idx) > 1:
            feat_mat[i, 0:6] = np.mean(np.abs(auc_features[:, idx]), axis=1).flatten()
            feat_mat[i, 6:19] = np.mean(pct_features[:, idx], axis=1).flatten()
            mod_mat[i, 0:13] = np.mean(mod_features[:, idx], axis=1).flatten()
        elif np.count_nonzero(idx) == 1:
            feat_mat[i, 0:6] = np.abs(auc_features[:, idx]).flatten()
            feat_mat[i, 6:19] = pct_features[:, idx].flatten()
            mod_mat[i, 0:13] = mod_features[:, idx].flatten()

        mod_mat[i, 13] = np.count_nonzero(idx)
    mod_mat[:, 13] = mod_mat[:, 13] / np.max(mod_mat[:, 13]) * 1.5
    # %% original scale ref

    su_ids = np.load(os.path.join(path, 'spike_clusters.npy'))
    spkTS = np.load(os.path.join(path, 'spike_times.npy'))

    su_ids = su_ids[np.squeeze((spkTS >= 30000 * 3600) & (spkTS < 30000 * (3600 + 1200)))]
    spkTS = spkTS[np.squeeze((spkTS >= 30000 * 3600) & (spkTS < 30000 * (3600 + 1200)))]

    id_set = np.unique(unitInfo['id'])
    # fh = plt.figure(figsize=(4.8, 9.6))
    fh, (ax1, ax2, ax3) = plt.subplots(1, 3)

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
        ax1.scatter(xs, ys, s=8, c='k', marker='|', edgecolors='none', alpha=0.002)
    ax1.set_ylim((0, 3840))
    ax1.set_xticks([])
    span = np.nanmax(feat_mat, axis=0) - np.nanmin(feat_mat, axis=0)
    span[span == 0] = 1
    feat_mat = (feat_mat - np.nanmin(feat_mat, axis=0)) / span
    im2 = ax2.imshow(feat_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=1)
    ax2.set_xticks([])
    ax2.set_yticks([])

    im3 = ax3.imshow(mod_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=1.5)
    [ax3.axvline(x, color='w', ls=':') for x in [2.5, 8.5]]
    ax3.set_xticks([])
    ax3.set_yticks([])
    fh.set_size_inches((4.8, 9.6))
    fh.savefig(os.path.join(path, 'scale.png'), dpi=300, bbox_inches="tight")
    plt.close(fh)


if __name__ == "__main__":
    DEBUGGING = False

    if DEBUGGING:
        path = r'K:\neupix\DataSum\191015-DPA-Learning2_29_g0\191015-DPA-Learning2_29_g0_imec0_cleaned'
        process_one_path(path)

    else:
        from multiprocessing import Pool

        curr_pool = Pool(processes=12)

        dpath = align.get_root_path()

        all_proc = []
        for path in zpy.traverse(dpath):
            print(path)
            all_proc.append(curr_pool.apply_async(process_one_path, args=(path,)))

        for one_proc in all_proc:
            one_proc.get()
            print('get')
