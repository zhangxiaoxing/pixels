# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:29:29 2019

@author: Libra
"""

import os
import selectivity as zpy
import phylib.utils._misc as phyutil

#%% Walk through

target='PL'
target_count=0
for path in zpy.traverse("K:/neupix/DataSum/"):
    
#%% filter for brain region
    FR_Th=1.0
    (bs_id,time_s,who)=zpy.get_bsid_duration_who(path)
    (mice_id, date, imec_no)=zpy.get_miceid_date_imecno(path)
  
    
    regionL=zpy.get_region_list();
    
    depthL=zpy.getTrackRegion(regionL, mice_id, date, imec_no, who)
    
    trials=zpy.get_trials(path)
    if trials is None:
        continue
    (perfType, perfIdx, selTrials) = zpy.judgePerformance(trials)
    
    spkNThresh = time_s * FR_Th
    unitInfo = phyutil.read_tsv(os.path.join(path, "cluster_info.tsv"))
    for row in unitInfo:
        if row["KSLabel"] == "good" and row["n_spikes"] >= spkNThresh:
            suDepth = row["depth"]
            fullR = zpy.matchDepth(suDepth, depthL)
            reg = zpy.combineSubRegion(fullR)
            if reg==target:
                # breakpoint()
                print(path)
                target_count+=1
        

#%% SPK->FR
                
                
#%% trial info -> filter trials




