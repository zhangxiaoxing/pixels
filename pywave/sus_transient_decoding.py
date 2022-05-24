"""
This code calculates cross time decoding and cross delay duration decoding, generates heatmaps
depends on the per_second_stats module for sustained and transient classification


"""

import sys
import per_second_stats
import os
import numpy as np

from multiprocessing import Pool
from ctd_actual import ctd_actual
from plot_ctd import plot_ctd

def cross_time_decoding(denovo=False, to_plot=False, delay=6, cpu=30, repeats=1000, wrs_thres=0, decoder='MCC'):
    if denovo:
        fstr = np.load("ctd.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        wrs_p = np.ones(len(features_per_su))
        sus_trans_flag = per_second_stats.process_all(denovo=True,
                                                      delay=delay)  # 33172 x 4, sust,trans,switch,unclassified
        sus_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[0, :], wrs_p > wrs_thres))[0]]
        trans_feat = [features_per_su[i] for i in
                      np.nonzero(np.logical_and(sus_trans_flag[1, :], wrs_p > wrs_thres))[0]]
        sus50 = []
        trans50 = []
        trans1000 = []
        sust_proc = []
        trans50_proc = []
        trans1000_proc = []
        if cpu > 1:
            curr_pool = Pool(processes=cpu)
            for i in range(np.ceil(repeats / 20).astype(np.int)):
                sust_proc.append(curr_pool.apply_async(ctd_actual, args=(sus_feat,),
                                                       kwds={"n_neuron": 50, "n_trial": (20, 25), "delay": delay,
                                                             "bin_range": np.arange(4, 52), "decoder": decoder}))

                trans50_proc.append(curr_pool.apply_async(ctd_actual, args=(trans_feat,),
                                                          kwds={"n_neuron": 50, "n_trial": (20, 25), "delay": delay,
                                                                "bin_range": np.arange(4, 52), "decoder": decoder}))
                trans1000_proc.append(curr_pool.apply_async(ctd_actual, args=(trans_feat,),
                                                            kwds={"n_neuron": 1000, "n_trial": (20, 25), "delay": delay,
                                                                  "bin_range": np.arange(4, 52), "decoder": decoder}))
            for one_proc in sust_proc:
                sus50.append(one_proc.get())

            for one_proc in trans50_proc:
                trans50.append(one_proc.get())

            for one_proc in trans1000_proc:
                trans1000.append(one_proc.get())

            curr_pool.close()
            curr_pool.join()
        else:
            sus50.append(ctd_actual(sus_feat, n_neuron=50, n_trial=(20, 25), delay=delay, bin_range=np.arange(4, 52),
                                    decoder=decoder))
            trans50.append(
                ctd_actual(trans_feat, n_neuron=50, n_trial=(20, 25), delay=delay, bin_range=np.arange(4, 52),
                           decoder=decoder))
            trans1000.append(
                ctd_actual(trans_feat, n_neuron=1000, n_trial=(20, 25), delay=delay, bin_range=np.arange(4, 52),
                           decoder=decoder))

        np.savez_compressed(f'sus_trans_ctd_{delay}_{repeats}.npz', sus50=sus50, trans50=trans50,
                            trans1000=trans1000)
    else:
        fstr = np.load(os.path.join('ctd', f'sus_trans_ctd_{delay}_1000.npz'), 'r')
        sus50 = fstr['sus50']
        trans50 = fstr['trans50']
        trans1000 = fstr['trans1000']

    if to_plot:
        plot_ctd(sus50,trans50,trans1000,delay,repeats)
        return (sus50, trans50, trans1000)



# %% main
if __name__ == "__main__":

    repeat = 10
    cpus = 1
    decoder = 'SVC'
    cross_time_decoding(denovo=True, to_plot=True, delay=6, cpu=cpus, repeats=repeat,
                        wrs_thres=0, decoder=decoder)
    cross_time_decoding(denovo=True, to_plot=True, delay=3, cpu=cpus, repeats=repeat,
                        wrs_thres=0, decoder=decoder)
    sys.exit()