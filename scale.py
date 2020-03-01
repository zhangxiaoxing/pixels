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


dpath = None
if os.path.exists("/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"):
    dpath = "/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"
elif os.path.exists(r"K:\neupix\DataSum"):
    dpath = r"K:\neupix\DataSum"
else:
    dpath = r"D:\neupix\DataSum"
for path in zpy.traverse(dpath):
    print(path)

    unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), sep="\t")

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
