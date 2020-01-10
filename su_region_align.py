# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:29:29 2019

@author: Libra
"""

import os
import h5py
import numpy as np
import pandas as pd
import selectivity as zpy


#%% Walk through

regionL=zpy.getRegionList();


for path in zpy.traverse("K:/neupix/DataSum/"):

    if os.path.isfile(os.path.join(path,'su_id2reg.csv')):
        continue
    if not os.path.isfile(os.path.join(path, "FR_All.hdf5")):
        print("missing one FR_All file, path: ", path)
        continue
    
    (bs_id,time_s,who)=zpy.get_bsid_duration_who(path)
    (mice_id, date, imec_no)=zpy.get_miceid_date_imecno(path)
    depthL=zpy.getTrackRegion(regionL, mice_id, date, imec_no, who)
    su_ids=[]
    su_region_corr=[]
    with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
        dset = ffr["SU_id"]
        su_ids = np.squeeze(np.array(dset, dtype="uint16"))
        
    unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), sep="\t")
    
    
    for one_su in su_ids:
        depth=unitInfo.loc[unitInfo['id']==one_su,['depth']].iat[0,0]
        reg=zpy.matchDepth(depth,depthL, date, mice_id, imec_no)
        su_region_corr.append([one_su,reg])
    
    tbl=pd.DataFrame(su_region_corr,columns=['index','region']).set_index('index')['region']
    tbl.to_csv(os.path.join(path,'su_id2reg.csv'),header=True)
        




