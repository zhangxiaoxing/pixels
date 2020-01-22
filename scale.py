# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 13:48:45 2020

@author: Libra
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


    
path=os.getcwd()
unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), sep="\t")

su_ids=np.load('spike_clusters.npy')
spkTS=np.load('spike_times.npy')

su_ids=su_ids[(spkTS>=30000*3600) & (spkTS<30000*(3600+2400))]
spkTS=spkTS[(spkTS>=30000*3600) & (spkTS<30000*(3600+2400))]

id_set=np.unique(unitInfo['id'])
fh=plt.figure(figsize=(4.8,9.6))
depth_list=[]
for one_id in id_set:
    depthDF=unitInfo.loc[unitInfo['id']==one_id,'depth']
    if depthDF.empty:
        print(one_id)
        continue
    depth=depthDF.iat[0]
    depth_list.append([one_id,depth])
    xs=spkTS[su_ids==one_id]/30000
    ys=np.ones_like(xs)*depth
    plt.scatter(xs,ys,s=2,c='k',marker='|',edgecolors='none',alpha=0.02)
fh.savefig(os.path.join(path, 'scale.png'), dpi=300, bbox_inches="tight")    