# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 23:34:48 2020

@author: Libra

Currently used only to plot 3s delay and 6s delay differential WM coding
"""



import os
import h5py
# import re
import scaling_stats
import numpy as np
import csv
import selectivity as zpy
import matplotlib.pyplot as plt


features_per_su=[]
reg_list=[]
for path in zpy.traverse("D:/neupix/DataSum/"):
    print(path)
    
    SU_ids = []
    trial_FR = []
    trials = []
    if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
        continue    
    done_read=False
    while not done_read:
        try:
            with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
                # print(list(ffr.keys()))
                if not "SU_id" in ffr.keys():
                    done_read=True
                    print("missing su_id key in path ", path)
                    continue
                # dset = ffr["SU_id"]
                # SU_ids = np.array(dset, dtype="uint16")
                dset = ffr["FR_All"]
                trial_FR = np.array(dset, dtype="double")
                dset = ffr["Trials"]
                trials = np.array(dset, dtype="double").T
            done_read=True
        except OSError:
            print("h5py read error handled")
        
    
    (perf_desc, perf_code, welltrain_window, correct_resp)=zpy.judgePerformance(trials)
    
    if perf_code!=3:
        continue
    
    suid_reg = []
    with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
        l = list(csv.reader(csvfile))[1:]
        suid_reg = [list(i) for i in zip(*l)]
        
        
    currStats=scaling_stats.scaling_stats()
    currStats.processGLMStats(trial_FR, trials, welltrain_window,correct_resp)
    features_per_su.append(currStats.get_features())
        
    reg_list.extend(suid_reg[1])
    

### save stats for future use
features=np.vstack(features_per_su)
reg_arr=np.array(reg_list)
np.savez_compressed('scaling_3s_6s.npz',features=features,reg_arr=reg_arr)    
    
### plot figure
fh=plt.figure(figsize=[4,3],dpi=300)
ax=fh.add_subplot(1,1,1)
h3=plt.plot(np.mean(features[:,:,0],axis=0),lw=1,c='r')
h6=plt.plot(np.mean(features[:,:,1],axis=0),lw=1,c='b')
[ax.axvline(x,c='k',ls=':',lw=0.5) for x in [12.5,16.5,28.5,32.5,40.5,44.5]]
ax.set_xticks([12.5, 32.5, 52.5])
ax.set_xticklabels([0,5,10])
ax.set_xlabel('Time (s)')
ax.set_ylabel('Sample selective at alpha=0.001')
ax.legend(['3s-delay','6s-delay'])    
plt.tight_layout(rect=[0, 0, 1, 0.95])
fh.savefig('scaling_3s_6s.png',dpi=300)
        
        
        
        
        
        