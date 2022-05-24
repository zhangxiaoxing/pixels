# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 11:48:13 2022

@author: Libra
"""

import h5py
import numpy as np
# from collections.abc import Iterable
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from proj_mat_hier import get_all_children, flatten
from get_allen_sets import get_conn_inj, get_conn_target

mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
structure_tree = mcc.get_structure_tree()
# grey_targets=list(flatten(get_all_children(structure_tree,567)))

grey_targets=list(flatten(get_all_children(structure_tree,567))) + list(flatten(get_all_children(structure_tree,343)))

inj_src=get_conn_inj()
inj_target=get_conn_target()

id_reg_map=structure_tree.get_id_acronym_map()
reg_id_map=structure_tree.get_name_map()

avail_inj=[x for x in grey_targets if x in inj_src]
avail_target=[x for x in grey_targets if x in inj_target]


src_mat=[]
src_reg_ids=[]
for src_id in avail_inj:
    print(src_id)
    one_src = structure_tree.get_structures_by_id([src_id])[0]
    one_src_exp = mcc.get_experiments(cre=False,injection_structure_ids=[one_src['id']])
    one_src_exp_ids=[exp['id'] for exp in one_src_exp]
    if not one_src_exp_ids:
        continue

    pm = mcc.get_projection_matrix(experiment_ids = one_src_exp_ids,
                                   projection_structure_ids = avail_target,
                                   hemisphere_ids= [3], # right hemisphere, ipsilateral
                                   parameter = 'projection_density')
    column_labels = [ c['label'] for c in pm['columns'] ]
    one_src_matrix = pm['matrix']
    src_mat.append(np.nanmedian(one_src_matrix,axis=0))
    src_reg_ids.append(src_id)

with h5py.File('proj_mat.hdf5','w') as fw:
    fw.create_dataset('grey_targets',data=np.array(avail_target).astype('int32'))
    fw.create_dataset('grey_srcs',data=np.array(src_reg_ids).astype('int32'))
    fw.create_dataset('src_target_matrix',data=np.array(src_mat).astype('float64'))