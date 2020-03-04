# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:38:49 2020

@author: Libra
"""
import os
import h5py
import csv
# import itertools
from GLM_stats import GLM_stats
from GLM_PCA_stats import GLM_PCA_stats
import numpy as np
from sklearn.decomposition import PCA
import su_region_align as align


def prepare_GLM():
    ### TODO: selective stats model
    curr_stats = GLM_PCA_stats()
    all_sess_list = []
    reg_list = []
    dpath = align.get_root_path()
    for path in align.traverse(dpath):
        print(path)
        # SU_ids = []
        trial_FR = None
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
                    # dset = ffr["SU_id"]
                    # SU_ids = np.array(dset, dtype="uint16")
                    dset = ffr["FR_All"]
                    trial_FR = np.array(dset, dtype="double")
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

        curr_stats.processGLMStats(trial_FR, trials, welltrain_window, correct_resp)

        onesession = curr_stats.getFeatures()
        all_sess_list.append(onesession)
        reg_list.extend(suid_reg[1])

    return (all_sess_list, reg_list)


# %% main
if __name__ == '__main__':
    (all_sess_list, reg_list) = prepare_GLM()
    ### save raw data file
    all_sess_arr = np.concatenate(tuple(all_sess_list), axis=1)
    reg_arr = np.array(reg_list)
    # np.savez_compressed('GLM_PCA_FR.npz',all_sess_arr=all_sess_arr,reg_arr=reg_arr)
    #
    #
    # ### load back saved raw data
    # fstr=np.load('GLM_PCA_FR.npz')
    # reg_arr=fstr['reg_arr']
    # all_sess_arr=fstr['all_sess_arr']

    pca = PCA(n_components=20)
    comp = pca.fit_transform(all_sess_arr)
    coeff = pca.components_
    ratio = pca.explained_variance_ratio_

    plot_figure = False

    reg_set = list(set(reg_arr.tolist()))
    reg_n = []
    for i in range(len(reg_set)):
        reg_n.append(f'{np.count_nonzero(reg_arr == reg_set[i])}, {reg_set[i]}')

    if plot_figure:
        import matplotlib.pyplot as plt

        fh = plt.figure(figsize=(10, 22), dpi=300)
    su_factors = []
    su_factors.append(reg_set)

    reg_count = []
    for reg in reg_set:
        reg_count.append(np.count_nonzero(reg_arr == reg))
    su_factors.append(reg_count)

    for sub in range(11):
        per_region_su_factor = np.zeros_like(reg_set, dtype=np.float64)
        for i in range(len(reg_set)):
            per_region_su_factor[i] = np.mean(all_sess_arr[sub, reg_arr == reg_set[i]])
        su_factors.append(per_region_su_factor.tolist())

        if plot_figure:
            feature_tag = ['sample selective during sample',
                           'sample selective during early delay',
                           'sample selective during late delay',
                           'sample selective ONLY during late delay',
                           'sample selective during decision making',
                           'test selective during decision making',
                           'pair selective during decision making',
                           'pair selective ONLY during reward',
                           'nonselective modulation during late delay',
                           'nonselective modulation during decision making',
                           'unmodulated']
            reg_idx = np.argsort(per_region_su_factor)
            ax = plt.subplot(11, 1, sub + 1)
            plt.bar(np.arange(len(reg_set)), per_region_su_factor[reg_idx])
            ax.set_xticks(np.arange(len(reg_set)))
            ax.set_xticklabels(np.array(reg_n)[reg_idx], rotation=90, ha='center', va='top', fontsize=7.5)
            ax.set_xlim(-1, len(reg_set))
            ax.set_ylabel(feature_tag[sub])

    if plot_figure:
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        fh.savefig("coding_feature.png", dpi=300, bbox_inches="tight")
        plt.show()
    su_factors = [list(i) for i in zip(*su_factors)]

    ### export csv for matlab GLM
    with open("PCA_Comp.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in su_factors:
            cwriter.writerow(row)
