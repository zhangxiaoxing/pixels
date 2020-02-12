# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:14:22 2020

@author: Libra
"""
import numpy as np
import scipy.stats as stats

class GLM_stats:
    def __init__(self):

        self.sample_sel_dur_sample=None
        self.sample_sel_dur_ED=None
        self.sample_sel_dur_LD=None
        self.sample_sel_only_LD=None
        self.sample_sel_dur_DM=None
        self.test_sel_dur_DM=None
        self.pair_sel_dur_DM=None
        self.pair_sel_only_reward=None
        self.non_sel_mod_dur_LD=None
        self.non_sel_mod_dur_DM=None
        self.non_mod=None   
        
        self.row_sel_6 = np.concatenate(
            (np.arange(16), np.arange(16, 40, 2), np.arange(40, 68))
        )
        self.row_sel_3 = np.arange(56)
        
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

    def processGLMStats(self, trial_FR, trials, welltrain_window=[],correct_resp=[]):

        ### TODO: when variables are empty
        FR_scale=np.concatenate((trial_FR[:,(trials[:,5]==3) & welltrain_window & correct_resp ,:][self.row_sel_3,:,:],
        trial_FR[:,(trials[:,5]==6) & welltrain_window & correct_resp,:][self.row_sel_6,:,:]),axis=1)
        
        trial_perf_sel=trials[((trials[:,5]==3) | (trials[:,5]==6)) & welltrain_window & correct_resp,:]
        
        trial_sel_left = trial_perf_sel[:, 2] == 4
 
        trial_sel_test_left = trial_perf_sel[:, 3] == 4
        
        pair_sel=(trial_perf_sel[:,2]!=trial_perf_sel[:,3])

        self.sample_sel_dur_sample=np.zeros((1,trial_FR.shape[2]))
        self.sample_sel_dur_ED=np.zeros_like(self.sample_sel_dur_sample)
        self.sample_sel_dur_LD=np.zeros_like(self.sample_sel_dur_sample)
        self.sample_sel_only_LD=np.zeros_like(self.sample_sel_dur_sample)
        self.sample_sel_dur_DM=np.zeros_like(self.sample_sel_dur_sample)
        self.test_sel_dur_DM=np.zeros_like(self.sample_sel_dur_sample)
        self.pair_sel_dur_DM=np.zeros_like(self.sample_sel_dur_sample)
        self.pair_sel_only_reward=np.zeros_like(self.sample_sel_dur_sample)
        self.non_sel_mod_dur_LD=np.zeros_like(self.sample_sel_dur_sample)
        self.non_sel_mod_dur_DM=np.zeros_like(self.sample_sel_dur_sample)
        self.non_mod=np.zeros_like(self.sample_sel_dur_sample)
        
        ### TODO: selective only during sample
        
        
        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(FR_scale[:, :, su_idx]).T
            # Uncomment if need baseline vector
            # (base_mean, base_std) = self.baselineVector(onesu)
            left_trials = onesu[trial_sel_left,:]
            right_trials = onesu[~trial_sel_left,:]

            left_test_trials = onesu[trial_sel_test_left,:]
            right_test_trials = onesu[~trial_sel_test_left,:]
            
            pair_trials=onesu[pair_sel,:]
            nonpair_trials=onesu[~pair_sel,:]
            
            # Sample 12:16, Delay 16:28, ED 18:22, LD 22:26, DM 28:36, reward 32:36
            
            self.sample_sel_dur_sample[0,su_idx]=self.bool_stats_test(left_trials[:,12:16],right_trials[:,12:16])
            self.sample_sel_dur_ED[0,su_idx]=self.bool_stats_test(left_trials[:,18:22],right_trials[:,18:22])
            self.sample_sel_dur_LD[0,su_idx]=self.bool_stats_test(left_trials[:,22:26],right_trials[:,22:26])
            self.sample_sel_only_LD[0,su_idx]=self.sample_sel_dur_LD[0,su_idx] and (not self.sample_sel_dur_sample[0,su_idx])
            self.sample_sel_dur_DM[0,su_idx]=self.bool_stats_test(left_trials[:,28:36],right_trials[:,28:36])
            
            self.test_sel_dur_DM[0,su_idx]=self.bool_stats_test(left_test_trials[:,28:36], right_test_trials[:,28:36])
            
            self.pair_sel_dur_DM[0,su_idx]=self.bool_stats_test(pair_trials[:,28:36], nonpair_trials[:,28:36])
            
            
            _pair_sel_dur_test=self.bool_stats_test(pair_trials[:,28:30], nonpair_trials[:,28:30])
            
            self.pair_sel_only_reward[0,su_idx]=(not _pair_sel_dur_test) and self.pair_sel_dur_DM[0,su_idx]
            
            self.non_sel_mod_dur_LD[0,su_idx]=((not self.sample_sel_dur_LD[0,su_idx]) and 
                                        self.bool_stats_test(onesu[:,4:8], onesu[:, 22:26]))

            self.non_sel_mod_dur_DM[0,su_idx]=((not self.pair_sel_dur_DM[0,su_idx]) and
                                        self.bool_stats_test(onesu[:,4:8], onesu[:, 30:34]))
            
            self.non_mod[0,su_idx]=((not self.bool_stats_test(onesu[:,4:8], onesu[:, 30:34])) &
                             (not self.bool_stats_test(onesu[:,4:8], onesu[:, 12:16])) &
                             (not self.bool_stats_test(onesu[:,4:8], onesu[:, 22:26])))
            
            ### TODO: selective only during sample

            
    def getFeatures(self):
        return np.concatenate((self.sample_sel_dur_sample,
            self.sample_sel_dur_ED,
            self.sample_sel_dur_LD,
            self.sample_sel_only_LD,
            self.sample_sel_dur_DM,
            self.test_sel_dur_DM,
            self.pair_sel_dur_DM,
            self.pair_sel_only_reward,
            self.non_sel_mod_dur_LD,
            self.non_sel_mod_dur_DM,
            self.non_mod
                ))