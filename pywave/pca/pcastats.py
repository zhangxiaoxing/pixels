# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:17:30 2021

@author: Libra
"""
import numpy as np
from sklearn.model_selection import KFold

class PCA_stats:
    def __init__(self):
        self.su_list = []
        self.avail = []
        self.kf = KFold(2, shuffle=True)

    def splitMean(self, FR, su_idx, mm, std, random_type=None):
        if random_type is None:
            out = np.squeeze(np.mean((FR[:, :, su_idx] - mm) / std, axis=1))
        elif random_type == 'half':
            sel = list(self.kf.split(np.arange(FR.shape[1])))[0]
            out = np.squeeze(np.mean((FR[:, sel[0], su_idx] - mm) / std, axis=1))
        else:
            sel = np.random.choice(FR.shape[1])
            out = np.squeeze((FR[:, sel, su_idx] - mm) / std)
        return out

        ### unprocessed data entry point

    def processCTDStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None, random_type=None):
        shifted_trial = np.vstack((np.zeros(6), trials[:-1, :]))

        ### TODO: when variables are empty

        FR_D3_S1_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 4, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D3_S2_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 8, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D3_S1_S2 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 4, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D3_S2_S2 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 8, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :, ]

        FR_D6_S1_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 4, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D6_S2_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 8, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]

        FR_D6_S1_S2 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 4, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D6_S2_S2 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 8, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :, ]

        for su_idx in range(trial_FR.shape[2]):
            (mm, std) = self.baselineVector(
                np.squeeze(trial_FR[:, np.logical_and(welltrain_window, correct_resp), su_idx])
            )
            if std > 0 and np.all([x.shape[1] >= 2 for x in (
                    FR_D3_S1_S1,
                    FR_D3_S1_S2,
                    FR_D3_S2_S1,
                    FR_D3_S2_S2,
                    FR_D6_S1_S1,
                    FR_D6_S1_S2,
                    FR_D6_S2_S1,
                    FR_D6_S2_S2,
            )]):
                self.su_list.append(
                    {
                        "S1_S1_3m": self.splitMean(FR_D3_S1_S1, su_idx, mm, std, random_type),
                        "S1_S2_3m": self.splitMean(FR_D3_S1_S2, su_idx, mm, std, random_type),
                        "S2_S1_3m": self.splitMean(FR_D3_S2_S1, su_idx, mm, std, random_type),
                        "S2_S2_3m": self.splitMean(FR_D3_S2_S2, su_idx, mm, std, random_type),
                        "S1_S1_6m": self.splitMean(FR_D6_S1_S1, su_idx, mm, std, random_type),
                        "S1_S2_6m": self.splitMean(FR_D6_S1_S2, su_idx, mm, std, random_type),
                        "S2_S1_6m": self.splitMean(FR_D6_S2_S1, su_idx, mm, std, random_type),
                        "S2_S2_6m": self.splitMean(FR_D6_S2_S2, su_idx, mm, std, random_type),
                    }
                )
                self.avail.append(True)
            else:
                self.avail.append(False)
    def baselineVector(self,one_su_trial_FR):
        base = one_su_trial_FR[2:10, :].flatten()
        if np.std(base):
            return (np.mean(base), np.std(base))
        else:
            base = one_su_trial_FR[-10:-1, :].flatten()
            if np.std(base):
                return (np.mean(base), np.std(base))
            else:
                print('flat base')
                base = one_su_trial_FR.flatten()
                return (0, np.std(base))

    def get_features(self):
        # breakpoint()
        return (self.su_list, self.avail)