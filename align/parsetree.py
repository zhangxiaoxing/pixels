# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 08:55:38 2021

@author: Libra
"""

import pandas as pd
import re
import numpy as np


def get_struct_tree():
    tree_file = r"K:\neupix\track_meta\structure_tree_safe_2017.csv"
    treetbl = pd.read_csv(tree_file,
                          usecols=['id','acronym','depth','structure_id_path',],
                          dtype={'id':'UInt32',
                                 'acronym':'string',
                                 'depth':'UInt8',
                                 'structure_id_path':'string'},
                          index_col='id')
    return treetbl

def get_tree_path(regstr):
    global treetbl

    regidx=treetbl.iloc[np.where(treetbl['acronym']==regstr)[0]].index[0]
    reg_depth=treetbl.loc[regidx,['depth']][0]
    reg_tree_path=treetbl.loc[regidx,['structure_id_path']][0]
    tree_idces=[int(x) for x in filter(None,reg_tree_path.split(r'/'))][2:]
    tree_str=[treetbl.loc[x,['acronym']][0] for x in tree_idces]
    return (regidx,reg_depth,tree_idces,tree_str)



if __name__=='parsetree':
    treetbl=get_struct_tree()