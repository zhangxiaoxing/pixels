# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 10:54:25 2021

@author: Libra
"""
import h5py
import numpy as np

def exporth5py(dict_stats):
    delay=dict_stats['delay']

    reg_depth=[int(dict_stats['reg'][x][1]) for x in range(len(dict_stats['reg']))]
    reg_tree=[dict_stats['reg'][x][2:] for x in range(len(dict_stats['reg']))]

    with h5py.File(f'transient_{delay}.hdf5', "w") as fw:
        fw.create_dataset('cluster_id', data=np.array(dict_stats['cid']).astype('uint16'))
        fw.create_dataset('selectivity', data=np.array(dict_stats['sel']).astype('float64'))
        fw.create_dataset('wrs_p',data=np.array(dict_stats['wrsp']).astype('float64'))
        fw.create_dataset('reg_tree_depth', data=np.array(reg_depth).astype('int8'))
        fw.create_dataset('reg_tree', data=np.array(reg_tree).astype('S10'))
        fw.create_dataset('path', data=np.array(dict_stats['folder']).astype('S120'))
        fw.create_dataset('delay', data=np.array(dict_stats['delay']).astype('uint8'))
        fw.create_dataset('trial_counts', data=np.array(dict_stats['trial_counts']).astype('uint8'))
        fw.create_dataset('wf_good', data=np.array(dict_stats['WF_good']).astype('int8'))