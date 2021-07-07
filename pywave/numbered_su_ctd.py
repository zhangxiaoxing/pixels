# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 09:52:59 2021

@author: Libra
"""

"""
This code calculates cross time decoding and cross delay duration decoding, generates heatmaps
depends on the per_second_stats module for sustained and transient classification


"""

import sys
import os
import numpy as np
import datautils

from multiprocessing import Pool
from ctd_actual import ctd_actual
from plot_ctd import plot_ctd

def numbered_su_ctd(denovo=False, to_plot=False, delay=6, cpu=30, repeats=1000, wrs_thres=0, decoder='MCC'):
    if denovo:
        #TODO retrive data from common origin
        # fstr = np.load("ctd.npz", allow_pickle=True)
        #TODO ramdom sample number of cells disregard of neuron type
        # features_per_su = fstr["features_per_su"].tolist()
        # wrs_p = np.ones(len(features_per_su))
        # sus_trans_flag = per_second_stats.process_all(denovo=True,
        #                                               delay=delay)  # 33172 x 4, sust,trans,switch,unclassified
        # sus_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[0, :], wrs_p > wrs_thres))[0]]
        # trans_feat = [features_per_su[i] for i in
        #               np.nonzero(np.logical_and(sus_trans_flag[1, :], wrs_p > wrs_thres))[0]]

        mixed_feat=datautils.get_dataset()

        su10 = []
        su100 = []
        su1000 = []
        sust_proc = []
        trans50_proc = []
        trans1000_proc = []
        if cpu > 1:
            curr_pool = Pool(processes=cpu)
            for i in range(np.ceil(repeats / 20).astype(np.int)):
                sust_proc.append(curr_pool.apply_async(ctd_actual, args=(mixed_feat,),
                                                       kwds={"n_neuron": 10, "n_trial": (20, 25), "delay": delay,
                                                             "bin_range": np.arange(4, 52), "decoder": decoder}))

                trans50_proc.append(curr_pool.apply_async(ctd_actual, args=(mixed_feat,),
                                                          kwds={"n_neuron": 100, "n_trial": (20, 25), "delay": delay,
                                                                "bin_range": np.arange(4, 52), "decoder": decoder}))
                trans1000_proc.append(curr_pool.apply_async(ctd_actual, args=(mixed_feat,),
                                                            kwds={"n_neuron": 1000, "n_trial": (20, 25), "delay": delay,
                                                                  "bin_range": np.arange(4, 52), "decoder": decoder}))
            for one_proc in sust_proc:
                su10.append(one_proc.get())

            for one_proc in trans50_proc:
                su100.append(one_proc.get())

            for one_proc in trans1000_proc:
                su1000.append(one_proc.get())

            curr_pool.close()
            curr_pool.join()
        else:
            su10.append(ctd_actual(mixed_feat, n_neuron=10, n_trial=(20, 25), delay=delay, bin_range=np.arange(4, 52),
                                    decoder=decoder))
            su100.append(
                ctd_actual(mixed_feat, n_neuron=100, n_trial=(20, 25), delay=delay, bin_range=np.arange(4, 52),
                           decoder=decoder))
            su1000.append(
                ctd_actual(mixed_feat, n_neuron=1000, n_trial=(20, 25), delay=delay, bin_range=np.arange(4, 52),
                           decoder=decoder))

        np.savez_compressed(f'sus_trans_ctd_{delay}_{repeats}.npz', su10=su10, su100=su100,su1000=su1000)
    else:
        fstr = np.load(os.path.join('ctd', f'sus_trans_ctd_{delay}_1000.npz'), 'r')
        su10 = fstr['su10']
        su100 = fstr['su100']
        su1000 = fstr['su1000']

    if to_plot:
        plot_ctd(su10,su100,su1000,delay,repeats)
        return (su10, su100, su1000)



# %% main
if __name__ == "__main__":

    repeat = 10
    cpus = 1
    decoder = 'SVC'
    numbered_su_ctd(denovo=True, to_plot=True, delay=6, cpu=cpus, repeats=repeat,
                        wrs_thres=0, decoder=decoder)
    numbered_su_ctd(denovo=True, to_plot=True, delay=3, cpu=cpus, repeats=repeat,
                        wrs_thres=0, decoder=decoder)
    sys.exit()
