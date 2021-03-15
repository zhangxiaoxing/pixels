# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:14:22 2020

@author: Libra


@author: Libra


"""


import platform
import sys
import pickle

if platform.system() == 'Windows':
    sys.path.insert(0, r"k:\code")

from per_sec_stats.prepare_data import prepare_data
from per_sec_stats.export import exporth5py

import align.su_region_align as align

def gen_align_files():
    align.gen_align_files()


if __name__ == "__main__":
    # prepare_data_sync()
    delay = 6
    if False:
        dict_stats=prepare_data(delay = delay, debug = True)
        pickle.dump(dict_stats,open(f'per_sec_sel_{delay}.p','wb'))
    else:
        dict_stats = pickle.load(open(f'per_sec_sel_{delay}.p','rb'))

    exporth5py(dict_stats)