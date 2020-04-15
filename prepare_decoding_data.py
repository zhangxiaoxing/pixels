# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:49:17 2020

@author: Libra
"""
import os

# import sys
import csv
import h5py
import numpy as np
import scipy.stats as stats
import su_region_align as align


class ctd_stats:
    def __init__(self):
        self.su_list = []
        self.use_ranksum = True

    def bool_stats_test(self, A, B):
        ### TODO alternative use of perm_test
        if self.use_ranksum:
            try:
                (_stat, p) = stats.mannwhitneyu(A.flatten(), B.flatten(), alternative="two-sided")
                return p < 0.001

            except ValueError:
                return False
        else:
            return self.exact_mc_perm_test(A, B, 1000)

    def exact_mc_perm_test(self, xs, ys, nmc):
        n, k = len(xs), 0
        diff = np.abs(np.mean(xs) - np.mean(ys))
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
        return k / nmc

    ### unprocessed data entry point

    def processCTDStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None):

        ### TODO: when variables are empty
        # [bin:trial:SU]
        # breakpoint()
        FR_D3_S1 = trial_FR[:,
                   np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 4, welltrain_window, correct_resp,)), axis=0, ),
                   :, ]
        FR_D3_S2 = trial_FR[:,
                   np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 8, welltrain_window, correct_resp,)), axis=0, ),
                   :, ]
        FR_D6_S1 = trial_FR[:,
                   np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 4, welltrain_window, correct_resp,)), axis=0, ),
                   :, ]
        FR_D6_S2 = trial_FR[:,
                   np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 8, welltrain_window, correct_resp,)), axis=0, ),
                   :, ]

        FR_D3_S1_ERR = trial_FR[:,
                       np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 4, np.logical_not(correct_resp))),
                              axis=0, ), :, ]
        FR_D3_S2_ERR = trial_FR[:,
                       np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 8, np.logical_not(correct_resp))),
                              axis=0, ), :, ]
        FR_D6_S1_ERR = trial_FR[:,
                       np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 4, np.logical_not(correct_resp))),
                              axis=0, ), :, ]
        FR_D6_S2_ERR = trial_FR[:,
                       np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 8, np.logical_not(correct_resp))),
                              axis=0, ), :, ]

        for su_idx in range(trial_FR.shape[2]):
            self.su_list.append(
                {
                    "S1_3": np.squeeze(FR_D3_S1[:, :, su_idx]),
                    "S2_3": np.squeeze(FR_D3_S2[:, :, su_idx]),
                    "S1_6": np.squeeze(FR_D6_S1[:, :, su_idx]),
                    "S2_6": np.squeeze(FR_D6_S2[:, :, su_idx]),
                    "S1_3_ERR": np.squeeze(FR_D3_S1_ERR[:, :, su_idx]),
                    "S2_3_ERR": np.squeeze(FR_D3_S2_ERR[:, :, su_idx]),
                    "S1_6_ERR": np.squeeze(FR_D6_S1_ERR[:, :, su_idx]),
                    "S2_6_ERR": np.squeeze(FR_D6_S2_ERR[:, :, su_idx]),
                }
            )

    def processPrevTrialStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None):

        ### TODO: when variables are empty
        # [bin:trial:SU]
        # breakpoint()
        
        shifted_trial=np.vstack((np.zeros(6),trials[1:,:]))
        
        FR_D3_S1 = trial_FR[:,
                   np.all(np.vstack((shifted_trial[:, 5] == 3, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)), axis=0, ),
                   :, ]
        FR_D3_S2 = trial_FR[:,
                   np.all(np.vstack((shifted_trial[:, 5] == 3, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)), axis=0, ),
                   :, ]
        FR_D6_S1 = trial_FR[:,
                   np.all(np.vstack((shifted_trial[:, 5] == 6, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)), axis=0, ),
                   :, ]
        FR_D6_S2 = trial_FR[:,
                   np.all(np.vstack((shifted_trial[:, 5] == 6, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)), axis=0, ),
                   :, ]

        FR_D3_S1_ERR = trial_FR[:,
                       np.all(np.vstack((shifted_trial[:, 5] == 3, shifted_trial[:, 2] == 4, np.logical_not(correct_resp))),
                              axis=0, ), :, ]
        FR_D3_S2_ERR = trial_FR[:,
                       np.all(np.vstack((shifted_trial[:, 5] == 3, shifted_trial[:, 2] == 8, np.logical_not(correct_resp))),
                              axis=0, ), :, ]
        FR_D6_S1_ERR = trial_FR[:,
                       np.all(np.vstack((shifted_trial[:, 5] == 6, shifted_trial[:, 2] == 4, np.logical_not(correct_resp))),
                              axis=0, ), :, ]
        FR_D6_S2_ERR = trial_FR[:,
                       np.all(np.vstack((shifted_trial[:, 5] == 6, shifted_trial[:, 2] == 8, np.logical_not(correct_resp))),
                              axis=0, ), :, ]

        for su_idx in range(trial_FR.shape[2]):
            self.su_list.append(
                {
                    "S1_3": np.squeeze(FR_D3_S1[:, :, su_idx]),
                    "S2_3": np.squeeze(FR_D3_S2[:, :, su_idx]),
                    "S1_6": np.squeeze(FR_D6_S1[:, :, su_idx]),
                    "S2_6": np.squeeze(FR_D6_S2[:, :, su_idx]),
                    "S1_3_ERR": np.squeeze(FR_D3_S1_ERR[:, :, su_idx]),
                    "S2_3_ERR": np.squeeze(FR_D3_S2_ERR[:, :, su_idx]),
                    "S1_6_ERR": np.squeeze(FR_D6_S1_ERR[:, :, su_idx]),
                    "S2_6_ERR": np.squeeze(FR_D6_S2_ERR[:, :, su_idx]),
                }
            )

    def get_features(self):
        # breakpoint()
        return self.su_list


### Main program entry point


def get_dataset(denovo=False, prevTrial=False):
    features_per_su = []
    reg_list = []
    if denovo:
        dpath = align.get_root_path()
        for path in align.traverse(dpath):
            print(path)
            trial_FR = None
            trials = None
            if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
                continue
            done_read = False
            while not done_read:
                try:
                    with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
                        if not "SU_id" in ffr.keys():
                            done_read = True
                            print("missing su_id key in path ", path)
                            break
                        dset = ffr["FR_All"]
                        trial_FR = np.array(dset, dtype="double")
                        dset = ffr["Trials"]
                        trials = np.array(dset, dtype="double").T
                    done_read = True
                except OSError:
                    print("h5py read error handled")

            if (trial_FR is None) or (trials is None):
                continue
            (_perf_desc, perf_code, welltrain_window, correct_resp,) = align.judgePerformance(trials, 75)

            if perf_code != 3:
                continue

            suid_reg = None
            with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
                l = list(csv.reader(csvfile))[1:]
                suid_reg = [list(i) for i in zip(*l)]

            currStats = ctd_stats()
            if prevTrial:
                currStats.processPrevTrialStats(trial_FR, trials, welltrain_window, correct_resp)
            else:
                currStats.processCTDStats(trial_FR, trials, welltrain_window, correct_resp)
            # as decoding is population non-linear statistic, will not calculate per neuron stats and average
            # will have to extract per neuron per trial stats and run population stats later
            # first get all neuron 25+25 trial
            # will need resample neurons later
            features_per_su.extend(currStats.get_features())
            reg_list.extend(suid_reg[1])

        ### DEBUG in small subsets
        # if len(features_per_su)>50:
        #     break

        ### save to npz file
        print("saving ctd data file...")
        if prevTrial:
            fn = 'prevTrial.npz'
        else:
            fn = 'ctd.npz'
        np.savez_compressed(fn, force_zip64=True, features_per_su=features_per_su, reg_list=reg_list)
        print("done saving")

### load from npz file
    else:
        fstr = np.load("ctd.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        reg_list = fstr["reg_list"].tolist()

    return (features_per_su, reg_list)


if __name__ == "__main__":
    (feat, regs) = get_dataset(denovo=True,prevTrial=True)
