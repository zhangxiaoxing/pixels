# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 14:03:42 2021

@author: Libra
"""
import h5py
import numpy as np
from collections.abc import Iterable
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache


def get_all_children(rootid,lvl=7):

    curr_str=structure_tree.get_structures_by_id([rootid])[0];
    if len(curr_str['structure_id_path'])<lvl:
        children=structure_tree.child_ids([rootid])[0]
        return [get_all_children(x) for x in children]
    else:
        return curr_str['id']


def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
structure_tree = mcc.get_structure_tree()
# FROM MOB
mob = structure_tree.get_structures_by_acronym(['MOB'])[0]
mob_exp = mcc.get_experiments(cre=False,injection_structure_ids=[mob['id']])
mob_exp_ids=[exp['id'] for exp in mob_exp]

ctx_targets=list(flatten(get_all_children(688)))

pm = mcc.get_projection_matrix(experiment_ids = mob_exp_ids,
                               projection_structure_ids = ctx_targets,
                               hemisphere_ids= [3], # right hemisphere, ipsilateral
                               parameter = 'projection_density')

column_labels = [ c['label'] for c in pm['columns'] ]
mob_matrix = pm['matrix']

#TO MO
mop = structure_tree.get_structures_by_acronym(['MOp'])[0]

mop_matrix=[]
mop_src=[]
for src in ctx_targets:
    print(src)
    mop_exp = mcc.get_experiments(cre=False,injection_structure_ids=[src])
    if not mop_exp:
        continue
    mop_exp_ids=[exp['id'] for exp in mop_exp]
    pm = mcc.get_projection_matrix(experiment_ids = mop_exp_ids,
                                   projection_structure_ids = [mop['id']],
                                   hemisphere_ids= [3], # right hemisphere, ipsilateral
                                   parameter = 'projection_density')
    mop_matrix.append(np.median(pm['matrix']))
    mop_src.append(src)

with h5py.File('proj_mat_hier.hdf5','w') as fw:
    fw.create_dataset('mob_targets',data=np.array(ctx_targets).astype('int32'))
    fw.create_dataset('mob_matrix',data=np.array(mob_matrix).astype('float64'))

    fw.create_dataset('mop_src',data=np.array(mop_src).astype('int32'))
    fw.create_dataset('mop_matrix',data=np.array(mop_matrix).astype('float64'))