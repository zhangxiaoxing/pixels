# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:57:04 2021

@author: Libra
"""

import wave.datautils as util
from pcastats import PCA_stats


def get_dataset():
    fr_per_su = util.get_dataset()
    currStats = PCA_stats()
    for onesu in fr_per_su:
        currStats.processCTDStats(
            onesu, random_type)
        (feat, avail) = currStats.get_features()
        features_per_su.extend(feat)
        avails.extend(avail)
        reg_list.extend(suid_reg[1])

    # DEBUG in small subsets
    # if len(features_per_su)>50:
    #     break

# save to npz file
    np.savez_compressed("ctd_prev_stats.npz",
                        features_per_su=features_per_su, reg_list=reg_list)

    return (features_per_su, reg_list, avails)