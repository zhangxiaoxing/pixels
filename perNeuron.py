
import os
import h5py
import re
import numpy as np
import selectivity as zpy
import pandas as pd
import matplotlib.pyplot as plt
from zxStats import zxStats
from scipy.cluster.hierarchy import dendrogram, ward



# def toFeature(trials,trial_FR,perf_judge):
#     S1T1Lick=perf_sel
    


su_hist=[]
last_path=[]
currStats = zxStats()
errorStats = zxStats()
currStats.regionName='All'
features_per_su=[]
error_features_per_su=[]
for path in zpy.traverse("K:/neupix/DataSum/"):
    
    SU_ids = []
    trial_FR = []
    trials = []
    with h5py.File(os.path.join(path, "FR_All.hdf5")) as ffr:
        # print(list(ffr.keys()))
        if not 'SU_id' in ffr.keys():
            print('missing su_id key in path ',path)
            continue
        dset = ffr["SU_id"]
        SU_ids = np.array(dset, dtype="uint16")
        dset = ffr["FR_All"]
        trial_FR = np.array(dset, dtype="double")
        dset = ffr["Trials"]
        trials = np.array(dset, dtype="double").T
        
    
    (perf_desc, perf_code, inWindow, correct_resp)=zpy.judgePerformance(trials)
    #  toReturn.extend(["wellTrained", 3, correctResp,welltrain_window])
    
    if perf_code!=3:
        continue
    
    currStats.addTrialFRs(trial_FR, trials, [], inWindow, correct_resp)
    errorStats.addTrialFRs(trial_FR, trials, [], inWindow, np.logical_not(correct_resp))
    
    
    features_per_su.append(currStats.getPerSUFeatureVector())    
    error_features_per_su.append(errorStats.getPerSUFeatureVector())    
    

    # if os.path.split(path)[0]==last_path:
    #     su_hist[-1]=su_hist[-1]+SU_ids.size
    # else:
    #     su_hist.append(SU_ids.size)
feat_per_su_arr = np.concatenate(tuple(features_per_su))
error_per_su_arr = np.concatenate(tuple(error_features_per_su))  

combinedFeatures=np.concatenate((feat_per_su_arr,error_per_su_arr),axis=1)
        
        
        
        
        
        
        