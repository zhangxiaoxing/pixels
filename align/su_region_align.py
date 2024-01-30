# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:29:29 2019

@author: Libra

This script is responsible for generating the su_id2reg.csv which is further used by other statistics


"""

import os, sys, re, glob, csv
import numpy as np
import pandas as pd
# sys.path.append('/home/zhangxx/code/pixels/align')
import fileutil as futil
import parsetree as parsetree

def imecNo2side(who_did, date, imecNo, mid): # assign probe # based on lab records
    if date == "20191130" and mid == "49":
        return "L"

    if date == "20191101" and mid == "26":
        return "R"

    if who_did == "HEM" and (int(date)) >= 20191028:
        if imecNo == "1":
            return "R"
        elif imecNo == "0":
            return "L"
    else:
        if imecNo == "1":
            return "L"
        elif imecNo == "0":
            return "R"

    print("Error parsing imec No")
    return "X"


def matchDepth(depth, depthL, date, mice_id, imec_no): # depth to brain-region-leaf
    if not depthL.empty:
        label = depthL.loc[
            (depthL["distance2tipLow"].to_numpy() <= depth) & (depth < depthL["distance2tipHigh"].to_numpy()),
            ["acronym"],
        ]
        if label.shape[0] == 1:
            return (label.iloc[0, 0],None)

    unlabeledRecord=[date, mice_id, imec_no, depth]
    return ('unlabeled',unlabeledRecord)


def getTrackRegionLN21(regionL, mice_id, date, imecNo): # all regions along a probe
    depthL = regionL.loc[
        (regionL["mouse_id"] == mice_id)
        & (regionL["implanting_date"] == date)
        & (regionL["imec"] == imecNo),
        ["acronym", "distance2tipLow", "distance2tipHigh"],
    ]
    return depthL

def getTrackRegion(regionL, mice_id, date, imecNo, who_did): # all regions along a probe
    depthL = regionL.loc[
        (regionL["mouse_id"] == mice_id)
        & (regionL["implanting_date"] == date)
        & (regionL["side"] == imecNo2side(who_did, date, imecNo, mice_id)),
        ["acronym", "distance2tipLow", "distance2tipHigh"],
    ]
    return depthL


def getRegionList(batch='WT'): # file interface
    if batch not in ('WT','LN21'):
        raise ValueError(f'Unknown experiment batch {batch}')
    if batch=='WT':
        site_file = r"/home/zhangxx/npdata/IHC00/NP_tracks_revisedNew.csv"
        regionL = pd.read_csv(site_file).astype(
            {"mouse_id": "str", "implanting_date": "str"}
        )[
            [
                "mouse_id",
                "implanting_date",
                "side",
                "acronym",
                "distance2tipHigh",
                "distance2tipLow",
            ]
        ]
    else:
        site_file = r"/home/zhangxx/npdata/IHC00/learning_tracks_new.CSV"
        regionL = pd.read_csv(site_file).astype(
                {"mouse_id": "str", "implanting_date": "str", "imec":"str"}
        )[
            [
                "mouse_id",
                "implanting_date",
                "imec",
                "acronym",
                "distance2tipHigh",
                "distance2tipLow",
            ]
        ]

    return regionL


# White paper coeff order
# ALM # Striatum # Thalamus # Midbrain # Medulla 

def reg2coeff(tree_str):
    # BS-HB-MY medulla
    # BS-HB-P midbrain
    # BS-IB-* thalamus
    # BS-MB-* midbrain
    # CB-*-* medulla
    # CH-CNU-*striatum
    # CH-CTX-* cortex
    if tree_str[1]=='CTX':
        return (0,'ALM')
    elif tree_str[1]=='CNU':
        return (1,'STR')
    elif tree_str[1]=='IB':
        return (2,'TH')
    elif tree_str[1]=='MB' or tree_str[2]=='P':
        return (3,'MB')
    elif tree_str[0]=='CB' or tree_str[2]=='MY':
        return (4,'MY')
    else:
        return (-1,'NONE')


def align_onefolder(unit_info,probe_path,regionL): # processing one session
        (mice_id, date, imec_no) = futil.get_miceid_date_imecno(probe_path) # meta data extracted from path and file tag
        if mice_id and date and imec_no:
            offset=int(imec_no)*10000
            if 'LN2' in probe_path:
                depthL = getTrackRegionLN21(regionL, mice_id, date, imec_no) # per probe brain regions
            else:
                ap_meta_path=glob.glob(probe_path[:-21]+"*.ap.meta")[0]
                (bs_id, time_s, who_did) = futil.get_bsid_duration_who(ap_meta_path) # base station id for lab records comparison
                depthL = getTrackRegion(regionL, mice_id, date, imec_no, who_did) # per probe brain regions

            su_region_corr = []

            dep_path=probe_path.replace('metrics.csv','channel_positions.npy')
            dep_lut=np.load(dep_path)
            map_path=probe_path.replace('metrics.csv','channel_map.npy')
            map_lut=np.load(map_path)
            unit_depth=np.array([dep_lut[np.where(map_lut==one_chan)[0],1] for one_chan in unit_info.peak_channel])

            unlabeled=[]
            ### sequentially match su ids with brain region tree
            for idx,one_su in unit_info.iterrows():
                depth = dep_lut[np.where(map_lut==one_su.peak_channel)[0],1]
                (reg,one_unlabeled) = matchDepth(depth, depthL, date, mice_id, imec_no)
                if reg=='unlabeled':
                    reg_depth=-1
                    tree_str=['',]
                else:
                    (regidx,reg_depth,tree_idces,tree_str)=parsetree.get_tree_path(reg)
                while len(tree_str)<6:
                    tree_str.append('')
                (coeffidx,coefftype)=reg2coeff(tree_str)
                su_region_corr.append([one_su.cluster_id+offset,coeffidx,coefftype,reg_depth]+tree_str[:6])
                if one_unlabeled:
                    unlabeled.append(one_unlabeled)
            return (su_region_corr,unlabeled)
        else:
#            print(f"error parsing data path \n{probe_path}")
            return([],[])

### traverse all folder
def gen_align_files(path):
#    unlabeledRecord = []
#    all_path=glob.glob(r'/home/zhangxx/npdata_out/**/*_imec?/',recursive=True)
#    sess_path=[re.findall(r'(^.*/)(.*?_imec[0-9])/',one_path)[0][0] for one_path in all_path]

#    for path in sess_path:
#         if os.path.exists(os.path.join(path, 'su_id2reg.csv')):
#             continue
#        if 'wyt' in path:
#            continue
    su_region_corr = []
    metrics_path=sorted(glob.glob('*/imec?_ks2/metrics.csv',root_dir=path))
    for one_probe in metrics_path:
        probe_path=os.path.join(path,one_probe)
        print(probe_path)
        metrics_df=pd.read_csv(probe_path)
        unit_info=metrics_df[["cluster_id","peak_channel"]]
        if 'LN2' in path:
            regionL = getRegionList('LN21')
        else:
            regionL = getRegionList('WT')

        (su_region,unlabeled)=align_onefolder(unit_info,probe_path,regionL)
        su_region_corr.extend(su_region)
#        unlabeledRecord.extend(unlabeled)

    if su_region_corr: # result export to file
        tbl = pd.DataFrame(su_region_corr,
                           columns=['index','coeff_idx','coeff_type','depth','d3','d4','d5','d6','d7','d8']).set_index('index')[:]
        tbl.to_csv(os.path.join(path, 'su_id2reg.csv'), header=True)

#    unlabeledRecord=list(filter(None, unlabeledRecord))
#    
#    with open("unlabeled.csv", "w", newline="") as f: # Incomplete data, developers only.
#        writer = csv.writer(f)
#        writer.writerows(unlabeledRecord)
#

if __name__=="__main__":
    if "snakemake" in globals():
        datapath=snakemake.wildcards[0]
        print(f"Processing path: {datapath}")
        gen_align_files(datapath)
