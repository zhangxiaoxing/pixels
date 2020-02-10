# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:29:29 2019

@author: Libra
"""

import os
import h5py
import numpy as np
import selectivity as zpy
import pandas as pd
import matplotlib.pyplot as plt
from zxStats import zxStats
from scipy.cluster.hierarchy import dendrogram, ward

# from sklearn.preprocessing import MinMaxScaler


#%% setup
plot_figures=False
# rois = ["PIR"]
rois = [
    "ACAd",
    "ACAv",
    "ACB",
    "AId",
    "AIv",
    "AON",
    "AUDd",
    "BST",
    "CA1",
    "CA3",
    "CEAl",
    "CEAm",
    "CLA",
    "CP",
    "DG-mo",
    "DP",
    "GPe",
    "HY",
    "IG",
    "ILA",
    "LA",
    "LH",
    "LP",
    "LSr",
    "MB",
    "MD",
    "MOp",
    "MOs",
    "OLF",
    "ORBl",
    "ORBm",
    "ORBvl",
    "PAL",
    "PIR",
    "PL",
    "PO",
    "PVT",
    "ProS",
    "RSPagl",
    "RSPd",
    "SI",
    "SSp-bfd",
    "SSp-n",
    "SSs",
    "STR",
    "TTd",
    "VISl",
    "VISp",
    "VPM",
    "aco",
    "ccb",
    "ccg",
    "ccs",
    "cing",
    "dhc",
    "fr",
    "or",
    "scwm",
]
error_rois = []
featured_rois = []
features = []
features_per_su=[]

for target in rois:
    # try:
    currStats = zxStats()
    currStats.regionName = target

    #%% loop through path
    debugCount = 0
    for path in zpy.traverse("K:/neupix/DataSum/"):

        debugCount += 1
        # if debugCount>50:
        #     break

        #%% filter for brain region
        if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
            print("Missing su_id to region matching file, path: ", path)
            continue

        su_id_regions = pd.read_csv(
            os.path.join(path, "su_id2reg.csv"), index_col="index"
        )

        su_id_list = su_id_regions[su_id_regions["region"] == target].index.to_numpy()

        if su_id_list.size == 0:
            continue

        SU_ids = []
        trial_FR = []
        trials = []
        with h5py.File(os.path.join(path, "FR_All.hdf5")) as ffr:
            # print(list(ffr.keys()))
            dset = ffr["SU_id"]
            SU_ids = np.array(dset, dtype="uint16")
            dset = ffr["FR_All"]
            trial_FR = np.array(dset, dtype="double")
            dset = ffr["Trials"]
            trials = np.array(dset, dtype="double").T

        (perf_desc, perf_code, inWindow, correct_resp) = zpy.judgePerformance(trials)
        if perf_code != 3:
            continue
        su_sel = np.squeeze(np.isin(SU_ids, su_id_list, assume_unique=True))

        currStats.addTrialFRs(trial_FR, trials, su_sel, inWindow, correct_resp)

    #%% plot actual figure
    if plot_figures:
        currStats.plotSummary()
    features.append(currStats.getFeatureVector())
    features_per_su.append(currStats.getPerSUFeatureVector())
    featured_rois.append(target)

# except:
#     error_rois.append(target)

feat_arr = np.array(features)
feat_per_su_arr = np.concatenate(tuple(features_per_su))
pd.Series(featured_rois, name="rois").to_csv("featured_rois.csv", header=True)
np.save("features.npy", feat_arr)
np.save("featuresPerSU.npy", feat_per_su_arr)

breakpoint()

# linear mapping to [0,1]

# 36, 20, 32 ,48, 48
segs = [range(36), range(36, 56), range(56, 88), range(88, 136), range(136, 192)]

feat_arr[:, 0:36] = (feat_arr[:, 0:36] - np.amin(feat_arr[:, 0:36])) / (
    np.amax(feat_arr[:, 0:36]) - np.amin(feat_arr[:, 0:36])
)
feat_arr[:, 36:56] = (feat_arr[:, 36:56] - np.amin(feat_arr[:, 36:56])) / (
    np.amax(feat_arr[:, 36:56]) - np.amin(feat_arr[:, 36:56])
)
feat_arr[:, 56:88] = (feat_arr[:, 56:88] - np.amin(feat_arr[:, 56:88])) / (
    np.amax(feat_arr[:, 56:88]) - np.amin(feat_arr[:, 56:88])
)
feat_arr[:, 88:136] = (feat_arr[:, 88:136] - np.amin(feat_arr[:, 88:136])) / (
    np.amax(feat_arr[:, 88:136]) - np.amin(feat_arr[:, 88:136])
)

feat_arr[:, 136:184] = (feat_arr[:, 136:184] - np.amin(feat_arr[:, 136:184])) / (
    np.amax(feat_arr[:, 136:184]) - np.amin(feat_arr[:, 136:184])
)

feat_arr[:, 184:] = (feat_arr[:, 184:] - np.amin(feat_arr[:, 184:])) / (
    np.amax(feat_arr[:, 184:]) - np.amin(feat_arr[:, 184:])
)


# scaler=MinMaxScaler(copy=False)
# scaler.fit_transform(feat_arr[:, 0:32])
# scaler.fit_transform(feat_arr[:, 32:64])
# scaler.fit_transform(feat_arr[:, 64:112])
# scaler.fit_transform(feat_arr[:, 112:])


#%% Feature Matrix
fh = plt.figure(figsize=(10, 7.5), dpi=300)
plt.imshow(feat_arr.T, cmap="jet", aspect="auto", vmin=0, vmax=1)
ax = plt.gca()
ax.set_xticks(range(len(featured_rois)))
ax.set_xticklabels(featured_rois, rotation="vertical")
ax.set_ylabel("Coding features")
plt.show()
fh.savefig("feature_map.png", dpi=300, bbox_inches="tight")

#%% Feature Matrix
fh = plt.figure(figsize=(100, 7.5), dpi=300)
plt.imshow(feat_per_su_arr.T, cmap="jet", vmin=0, vmax=1)
ax = plt.gca()
# ax.set_xticks(range(len(featured_rois)))
# ax.set_xticklabels(featured_rois, rotation="vertical")
ax.set_xlabel('# of single units')
ax.set_ylabel("Coding features")
plt.show()
fh.savefig("feature_map_per_su_100.png", dpi=300, bbox_inches="tight")




from opgene_tsne_map import opgene_tsne
myTsne=opgene_tsne()
myTsne.rois=rois
#36, 20, 32 ,48, 48, 48
feat_tsne=myTsne.calc_tsne(feat_arr[:,:136])
myTsne.plot_tsne_agglo(5)

myTsne.plot_feature_decomp()

myTsne.plot_opgene_mapping('Simple_DPA_DPA-DM-Hit',True)
myTsne.plot_opgene_mapping('Simple_DPA_DPA-DM-CR',True)
myTsne.plot_opgene_mapping('Simple_DPA_DPA-DM',True)

myTsne.plot_opgene_mapping('Simple_DPA_DPA-LD-Hit',True)
myTsne.plot_opgene_mapping('Simple_DPA_DPA-LD-CR',True)
myTsne.plot_opgene_mapping('Simple_DPA_DPA-LD',True)

myTsne.plot_opgene_mapping('Simple_DPA_DPA-ED-Hit',True)
myTsne.plot_opgene_mapping('Simple_DPA_DPA-ED-CR',True)
myTsne.plot_opgene_mapping('Simple_DPA_DPA-ED',True)



#%% dendrograph

linkage_array = ward(feat_arr)
fh = plt.figure(figsize=[6, 6], dpi=300)
dendrogram(linkage_array, labels=rois, color_threshold=3.5)
fh.savefig("dendrogram.png", dpi=300, bbox_inches="tight")



