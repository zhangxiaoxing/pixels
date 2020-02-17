# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:49:17 2020

@author: Libra
"""
import os
import h5py
import scaling_stats
import csv
import pickle
import numpy as np
import selectivity as zpy
import matplotlib.pyplot as plt
import scipy.stats as stats

class ctd_stats:
    def __init__(self):
        self.su_list=[]
        
        self.use_ranksum=True
        
    def bool_stats_test(self, A,B):
        ### TODO alternative use of perm_test
        if self.use_ranksum:
            try:
                (_stat, p) = stats.mannwhitneyu(
                    A.flatten(),
                    B.flatten(),
                    alternative="two-sided",
                )
                return p < 0.001
            
            except ValueError:
                return False
        else:
            return self.exact_mc_perm_test(A,B,1000)
                
            

    def exact_mc_perm_test(self, xs, ys, nmc):
        n, k = len(xs), 0
        diff = np.abs(np.mean(xs) - np.mean(ys))
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
        return k / nmc
    
    ### unprocessed data entry point

    def processGLMStats(self, trial_FR, trials, welltrain_window=None,correct_resp=None):

        ### TODO: when variables are empty
        # [bin:trial:SU]
        # breakpoint()
        FR_D3_S1=trial_FR[:,np.all(np.vstack((trials[:,5]==3, trials[:,2]==4 , welltrain_window , correct_resp)),axis=0) ,:]
        FR_D3_S2=trial_FR[:,np.all(np.vstack((trials[:,5]==3, trials[:,2]==8 , welltrain_window , correct_resp)),axis=0) ,:]
        FR_D6_S1=trial_FR[:,np.all(np.vstack((trials[:,5]==6, trials[:,2]==4 , welltrain_window , correct_resp)),axis=0) ,:]
        FR_D6_S2=trial_FR[:,np.all(np.vstack((trials[:,5]==6, trials[:,2]==8 , welltrain_window , correct_resp)),axis=0) ,:]
        
        FR_D3_S1_ERR=trial_FR[:,np.all(np.vstack((trials[:,5]==3, trials[:,2]==4 , np.logical_not(correct_resp))),axis=0) ,:]
        FR_D3_S2_ERR=trial_FR[:,np.all(np.vstack((trials[:,5]==3, trials[:,2]==8 , np.logical_not(correct_resp))),axis=0) ,:]
        FR_D6_S1_ERR=trial_FR[:,np.all(np.vstack((trials[:,5]==6, trials[:,2]==4 , np.logical_not(correct_resp))),axis=0) ,:]
        FR_D6_S2_ERR=trial_FR[:,np.all(np.vstack((trials[:,5]==6, trials[:,2]==8 , np.logical_not(correct_resp))),axis=0) ,:]

        for su_idx in range(trial_FR.shape[2]):
            self.su_list.append({'S1_3':np.squeeze(FR_D3_S1[:,:,su_idx]),
                                 'S2_3':np.squeeze(FR_D3_S2[:,:,su_idx]),
                                 'S1_6':np.squeeze(FR_D6_S1[:,:,su_idx]),
                                 'S2_6':np.squeeze(FR_D6_S2[:,:,su_idx]),
                                 'S1_3_ERR':np.squeeze(FR_D3_S1_ERR[:,:,su_idx]),
                                 'S2_3_ERR':np.squeeze(FR_D3_S2_ERR[:,:,su_idx]),
                                 'S1_6_ERR':np.squeeze(FR_D6_S1_ERR[:,:,su_idx]),
                                 'S2_6_ERR':np.squeeze(FR_D6_S2_ERR[:,:,su_idx]),
                                 })
    def get_features(self):
        # breakpoint()
        return self.su_list
    

### program entry point    
        
if __name__=='__main__':
    
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
            
        currStats=ctd_stats()
        currStats.processGLMStats(trial_FR, trials, welltrain_window,correct_resp)
        features_per_su.extend(currStats.get_features())
        reg_list.extend(suid_reg[1])
        
        ### test in small subsets
        # if len(features_per_su)>50:
        #     break

### save file to pickle    
with open('ctd.pickle','wb') as fh:
    pickle.dump({'features_per_su':features_per_su, 'reg_list':reg_list},fh)


### load back saved file on demand    
# dstr=None    
# with open('ctd.pickle','rb') as fh:
#     dstr=pickle.load(fh)