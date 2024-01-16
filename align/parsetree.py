# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 08:55:38 2021

@author: Libra
"""

import pandas as pd
import numpy as np


def _get_struct_tree(): #import from Allen CCF v3 file
    tree_file = r"/home/zhangxx/npdata/IHC00/structure_tree_safe.csv"
    # pandas #44079 work around
    return pd.read_csv(tree_file,
                          usecols=['id','acronym','depth','structure_id_path',],
                          dtype={'id':'UInt32',
                                  'acronym':'string',
                                  'depth':'UInt8',
                                  'structure_id_path':'string'},
                          ).set_index(['id',])



def get_tree_path(regstr): #build up entire region tree from leaf node
    try:
        regidx=_treetbl.iloc[np.where(_treetbl['acronym']==regstr)[0]].index[0]
    except Exception:
        breakpoint()
    reg_depth=_treetbl.loc[regidx]['depth']
    reg_tree_path=_treetbl.loc[regidx]['structure_id_path']
    tree_idces=[int(x) for x in filter(None,reg_tree_path.split(r'/'))][2:]
    tree_str=[_treetbl.loc[x]['acronym'] for x in tree_idces]
    return (regidx,reg_depth,tree_idces,tree_str)


# eval during import
_treetbl=_get_struct_tree()
