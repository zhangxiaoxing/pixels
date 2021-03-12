# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 08:45:24 2021

@author: Libra
"""


def process_select_stats(trial_FR, trials, welltrain_window=None, correct_resp=None, delay=6):
    # trial_FR [bin, trial, SU]
    ### TODO: when variables are none

    firing_rate_6 = trial_FR[:, (trials[:, 5] == 6) & welltrain_window & correct_resp, :]
    firing_rate_3 = trial_FR[:, (trials[:, 5] == 3) & welltrain_window & correct_resp, :]
    trial_perf_sel_6 = trials[(trials[:, 5] == 6) & welltrain_window & correct_resp, :]
    trial_perf_sel_3 = trials[(trials[:, 5] == 3) & welltrain_window & correct_resp, :]

    trial_sel_left_6 = trial_perf_sel_6[:, 2] == 4
    trial_sel_left_3 = trial_perf_sel_3[:, 2] == 4

    self.per_sec_selectvity = np.zeros((delay+1, trial_FR.shape[2]))
    self.per_sec_auc = np.zeros((delay+1, trial_FR.shape[2]))
    self.per_sec_fr = np.zeros((delay+1, trial_FR.shape[2]))
    self.per_sec_wrs_p = np.zeros((delay + 1, trial_FR.shape[2]))

    for su_idx in range(trial_FR.shape[2]):
        onesu_6 = np.squeeze(firing_rate_6[:, :, su_idx]).T
        onesu_3 = np.squeeze(firing_rate_3[:, :, su_idx]).T
        left_trials_6 = onesu_6[trial_sel_left_6, :]
        right_trials_6 = onesu_6[~trial_sel_left_6, :]

        left_trials_3 = onesu_3[trial_sel_left_3, :]
        right_trials_3 = onesu_3[~trial_sel_left_3, :]

        baseL = np.mean(trial_FR[:10, :, su_idx][:, (trials[:, 2] == 4) & welltrain_window & correct_resp],
                        axis=0)
        baseR = np.mean(trial_FR[:10, :, su_idx][:, (trials[:, 2] == 8) & welltrain_window & correct_resp],
                        axis=0)
        self.baseline_sel[:,su_idx] = self.bool_stats_test(baseL, baseR)


        for bin_idx in range(0, delay+1):
            bins = np.arange(bin_idx * 4 + 12, bin_idx * 4 + 16)
            if delay == 6:
                spksum = (np.sum(left_trials_6[:, bins]) + np.sum(right_trials_6[:, bins]))
                bsstat=self.bool_stats_test(onesu_6[:, bins], onesu_6[:, 6:10], delay+1)[0]
            else:
                spksum = (np.sum(left_trials_3[:, bins]) + np.sum(right_trials_3[:, bins]))
                bsstat=self.bool_stats_test(onesu_3[:, bins], onesu_3[:, 6:10], delay+1)[0]

            if spksum == 0:
                self.per_sec_sel_raw[bin_idx, su_idx] = 0
                self.per_sec_auc[bin_idx, su_idx] = 0.5
                self.per_sec_fr[bin_idx, su_idx] = 0
                self.per_sec_sel[bin_idx, su_idx] = False
                self.per_sec_wrs_p[bin_idx, su_idx] = 1
            else:
                if delay == 6:
                    binstat=self.bool_stats_test(left_trials_6[:, bins],right_trials_6[:, bins], delay+1)
                    self.per_sec_fr[bin_idx, su_idx] =np.mean(np.concatenate((left_trials_6[:,bins],right_trials_6[:,bins])))
                else:
                    binstat=self.bool_stats_test(left_trials_3[:, bins],right_trials_3[:, bins], delay+1)
                    self.per_sec_fr[bin_idx, su_idx] =np.mean(np.concatenate((left_trials_3[:,bins],right_trials_3[:,bins])))

                self.per_sec_sel_raw[bin_idx, su_idx] = binstat[2]
                self.per_sec_sel[bin_idx, su_idx] = binstat[0]
                self.per_sec_auc[bin_idx, su_idx] =binstat[3]
                self.per_sec_wrs_p[bin_idx, su_idx] = binstat[1]

            if self.per_sec_sel[bin_idx, su_idx]:
                if (delay == 6 and np.mean(left_trials_6[:, bins]) > np.mean(right_trials_6[:, bins]))\
                    or (delay == 3 and np.mean(left_trials_3[:, bins]) > np.mean(right_trials_3[:, bins])):

                        self.per_sec_prefS1[bin_idx, su_idx] = 1
                else:
                    self.per_sec_prefS2[bin_idx, su_idx] = 1



            self.non_sel_mod[bin_idx, su_idx] = bsstat and not self.per_sec_sel[bin_idx, su_idx]