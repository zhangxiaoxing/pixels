# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 15:48:09 2020

@author: Libra
"""


import os
import h5py
import numpy as np
import selectivity as zpy
import tensortools as tt
import matplotlib.pyplot as plt


def rearrange_row(trials, trial_FR):
    row_sel_6 = np.concatenate((np.arange(16), np.arange(16, 40, 2), np.arange(40, 68)))
    row_sel_3 = np.arange(56)
    row_arr_out = []
    for i in range(trials.shape[0]):
        if trials[i, 5] <= 3:
            row_arr_out.append(trial_FR[row_sel_3, i, :])
        elif trials[i, 5] >= 6:
            row_arr_out.append(trial_FR[row_sel_6, i, :])
        else:
            print("Error delay duration")

    return np.moveaxis(np.array(row_arr_out), 2, 0)


def rearrange_block(trials, upper_bound=0):
    if upper_bound == 0:
        upper_bound = trials.shape[0]
    block_begin = 0
    all_out = np.zeros(upper_bound, dtype=np.int32)
    while block_begin < upper_bound:
        unregistered = np.full(8, True)
        for i in range(block_begin, block_begin + 8):
            if trials[i, 2] == 4 and trials[i, 3] == 4:
                if unregistered[0]:
                    all_out[block_begin] = i
                    unregistered[0] = False
                else:
                    all_out[block_begin + 1] = i
                    unregistered[4] = False
            elif trials[i, 2] == 4 and trials[i, 3] == 8:
                if unregistered[1]:
                    all_out[block_begin + 2] = i
                    unregistered[1] = False
                else:
                    all_out[block_begin + 3] = i
                    unregistered[5] = False
            elif trials[i, 2] == 8 and trials[i, 3] == 4:
                if unregistered[2]:
                    all_out[block_begin + 4] = i
                    unregistered[2] = False
                else:
                    all_out[block_begin + 5] = i
                    unregistered[6] = False
            else:  # 8, 8
                if unregistered[3]:
                    all_out[block_begin + 6] = i
                    unregistered[3] = False
                else:
                    all_out[block_begin + 7] = i
                    unregistered[7] = False
        if np.any(unregistered):
            return (False, all_out)

        block_begin += 8

    return (True, all_out)


def rearrange_sep_block(trials, upper_bound=0):
    if upper_bound == 0:
        upper_bound = trials.shape[0]
    block_begin = 0
    S1Pair = []
    S1nonPair = []
    S2Pair = []
    S2nonPair = []
    while block_begin < upper_bound:
        unregistered = np.full(8, True)
        for i in range(block_begin, block_begin + 8):
            if trials[i, 2] == 4 and trials[i, 3] == 8:
                if unregistered[0]:
                    S1Pair.append(i)
                    unregistered[0] = False
                else:
                    S1Pair.append(i)
                    unregistered[4] = False
            elif trials[i, 2] == 4 and trials[i, 3] == 4:
                if unregistered[1]:
                    S1nonPair.append(i)
                    unregistered[1] = False
                else:
                    S1nonPair.append(i)
                    unregistered[5] = False
            elif trials[i, 2] == 8 and trials[i, 3] == 4:
                if unregistered[2]:
                    S2Pair.append(i)
                    unregistered[2] = False
                else:
                    S2Pair.append(i)
                    unregistered[6] = False
            else:  # 8, 8
                if unregistered[3]:
                    S2nonPair.append(i)
                    unregistered[3] = False
                else:
                    S2nonPair.append(i)
                    unregistered[7] = False
        if np.any(unregistered):
            return (False, [])

        block_begin += 8
    all_out = np.array(S1Pair + S1nonPair + S2Pair + S2nonPair)
    return (True, all_out)


def normalize(all_sess_arr):  # SU,trials, bins
    arrmin = np.amin(all_sess_arr, axis=(1, 2))
    arrspan = np.amax(all_sess_arr, axis=(1, 2)) - arrmin
    all_sess_arr = np.moveaxis(all_sess_arr, 0, 2)
    all_sess_arr -= arrmin
    all_sess_arr /= arrspan
    return np.moveaxis(all_sess_arr, 2, 0)


def nonneg_tca(X, R, prefix='',max_iter=500):

    # Fit CP tensor decomposition (two times).
    U = tt.ncp_bcd(X, rank=R, verbose=True, max_iter=max_iter,tol=1E-6)
    V = tt.ncp_bcd(X, rank=R, verbose=True, max_iter=max_iter,tol=1E-6)
    
    # Compare the low-dimensional factors from the two fits.
    # fig, ax, po = tt.plot_factors(U.factors)
    # tt.plot_factors(V.factors, fig=fig)
    # fig.suptitle("raw models")
    # fig.tight_layout()

    # Align the two fits and print a similarity score.
    sim = tt.kruskal_align(U.factors, V.factors, permute_U=True, permute_V=True)
    print(sim)

    
    
    
    # Plot the results again to see alignment.
    fig, ax, po = tt.plot_factors(U.factors, plots=["scatter", "scatter", "line"])
    tt.plot_factors(V.factors, plots=["scatter", "scatter", "line"], fig=fig)
    [x.set_xticks([11.5, 15.5, 27.5, 31.5]) for x in ax[:, 2]]

    ax[-1, 0].set_xlabel("SU #")
    ax[-1, 1].set_xlabel("Trial #")
    ax[-1, 2].set_xlabel("Time (s)")
    ax[-1, 2].set_xticklabels(["S", "+1", "T", "+1"])

    fig.suptitle("aligned models")
    fig.tight_layout()

    # Show plots.
    plt.show()
    fig.set_size_inches(40, 40)
    fig.set_dpi(300)
    fig.savefig(
        prefix+"nonneg_TCA_trial_" + str(X.shape[1]) + "_R" + str(R) + ".png",
        dpi=300,
        bbox_inches="tight",
    )
    return (U, V, sim)



def run_tca(trial_target, sep_blocks=False):
    ntrialsCount = []
    all_sess_list = []

    for path in zpy.traverse("K:/neupix/DataSum/"):
        print(path)
        # SU_ids = []
        trial_FR = []
        trials = []
        with h5py.File(os.path.join(path, "FR_All.hdf5")) as ffr:
            # print(list(ffr.keys()))
            if not "SU_id" in ffr.keys():
                print("missing su_id key in path ", path)
                continue
            # dset = ffr["SU_id"]
            # SU_ids = np.array(dset, dtype="uint16")
            dset = ffr["FR_All"]
            trial_FR = np.array(dset, dtype="double")
            dset = ffr["Trials"]
            trials = np.array(dset, dtype="double").T

        (perf_desc, perf_code, inWindow, correct_resp) = zpy.judgePerformance(trials)
        #  toReturn.extend(["wellTrained", 3, correctResp,welltrain_window])

        if perf_code != 3:
            continue
        ntrialsCount.append(trials.shape[0])

        if trials.shape[0] < trial_target:
            continue
        
        if sep_blocks:
            (reg_all, matched_index) = rearrange_sep_block(trials, trial_target)
            
        else:
            (reg_all, matched_index) = rearrange_block(trials, trial_target)

        if reg_all:
            merged = rearrange_row(trials, trial_FR)
            # onesession=merged[:,matched_index,:].reshape((merged.shape[0],-1))
            onesession = merged[:, matched_index, :]
            all_sess_list.append(onesession)

    all_sess_arr = np.concatenate(tuple(all_sess_list), axis=0)
    all_sess_arr = normalize(all_sess_arr)
    opti_param = []
    
    sep_str='sepblock_' if sep_blocks else 'consec_'
    for R in range(5, 16):
        (U, V, sim) = nonneg_tca(all_sess_arr, R, prefix=sep_str,max_iter=1000)
        np.savez_compressed(sep_blocks+"tensor_comp_trial" + str(all_sess_arr.shape[1]) + "_R" + str(R) + ".npz",
                            SU=U.factors.factors[0],
                            cross_trial=U.factors.factors[1],
                            in_trial=U.factors.factors[2])
        opti_param.append([R, U.obj, V.obj, sim])
    np.save(
        sep_str+"nonneg_trials" + str(trial_target) + "_opti_params.npy", np.array(opti_param)
    )
