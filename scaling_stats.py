# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 23:34:48 2020

@author: Libra
"""

import numpy as np
import scipy.stats as stats

class scaling_stats:
    def __init__(self):
        self.sample_sel_D3=None
        self.sample_sel_D6=None
        
        self.use_ranksum=True
        
    def bool_stats_test(self, A,B):
        ### TODO alternative use of perm_test
        if self.use_ranksum:
            try:
                (stat, p) = stats.mannwhitneyu(
                    A.flatten(),
                    B.flatten(),
                    alternative="two-sided",
                )
                return p < 0.001
            
            except ValueError:
                return False
            

    def exact_mc_perm_test(self, xs, ys, nmc):
        n, k = len(xs), 0
        diff = np.abs(np.mean(xs) - np.mean(ys))
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
        return k / nmc
    
    ### unprocessed data entry point

    def processGLMStats(self, trial_FR, trials, welltrain_window=[],correct_resp=[]):

        ### TODO: when variables are empty
        # [bin:trial:SU]
        # breakpoint()
        FR_D3_S1=trial_FR[:,np.all(np.vstack((trials[:,5]==3, trials[:,2]==4 , welltrain_window , correct_resp)),axis=0) ,:];
        FR_D3_S2=trial_FR[:,np.all(np.vstack((trials[:,5]==3, trials[:,2]==8 , welltrain_window , correct_resp)),axis=0) ,:];
        FR_D6_S1=trial_FR[:,np.all(np.vstack((trials[:,5]==6, trials[:,2]==4 , welltrain_window , correct_resp)),axis=0) ,:];
        FR_D6_S2=trial_FR[:,np.all(np.vstack((trials[:,5]==6, trials[:,2]==8 , welltrain_window , correct_resp)),axis=0) ,:];

        self.sample_sel_D3=np.zeros((trial_FR.shape[2],trial_FR.shape[0]))
        self.sample_sel_D6=np.zeros_like(self.sample_sel_D3)
        
        for su_idx in range(trial_FR.shape[2]):
            for bin_idx in range(trial_FR.shape[0]):
                # Uncomment if need baseline vector
                # (base_mean, base_std) = self.baselineVector(onesu)
                # if su_idx>=FR_D3_S1.shape[2]:
                    # breakpoint()
                self.sample_sel_D3[su_idx,bin_idx]=self.bool_stats_test(FR_D3_S1[bin_idx,:,su_idx],FR_D3_S2[bin_idx,:,su_idx])
                self.sample_sel_D6[su_idx,bin_idx]=self.bool_stats_test(FR_D6_S1[bin_idx,:,su_idx],FR_D6_S2[bin_idx,:,su_idx])
            
    def get_features(self):
        # breakpoint()
        return np.dstack((self.sample_sel_D3,
            self.sample_sel_D6
                ))