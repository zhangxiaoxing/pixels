# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 21:08:04 2020

@author: Libra
"""

import numpy as np
from pixelStats import baselineVector


class GLM_PCA_stats:
    def __init__(self):

        self.S1T1=None
        self.S1T2=None
        self.S2T1=None
        self.S2T2=None
        
        self.row_sel_6 = np.concatenate(
            (np.arange(16), np.arange(16, 40, 2), np.arange(40, 68))
        )
        self.row_sel_3 = np.arange(56)
        

    def processGLMStats(self, trial_FR, trials, welltrain_window=[],correct_resp=[]):

        ### TODO: when variables are empty
        
        FR_scale=np.concatenate((trial_FR[:,(trials[:,5]==3) & welltrain_window & correct_resp ,:][self.row_sel_3,:,:],
        trial_FR[:,(trials[:,5]==6) & welltrain_window & correct_resp,:][self.row_sel_6,:,:]),axis=1)
        
        trial_perf_sel=np.concatenate((trials[(trials[:,5]==3) & welltrain_window & correct_resp,:],
                                       trials[(trials[:,5]==6) & welltrain_window & correct_resp,:]),axis=0)
        
        # S1=4 S2=8 T1=8 T2 =4
        trial_S1T1 = (trial_perf_sel[:, 2] == 4) & (trial_perf_sel[:, 3] == 8)
        trial_S1T2 = (trial_perf_sel[:, 2] == 4) & (trial_perf_sel[:, 3] == 4)
        trial_S2T1 = (trial_perf_sel[:, 2] == 8) & (trial_perf_sel[:, 3] == 8)
        trial_S2T2 = (trial_perf_sel[:, 2] == 8) & (trial_perf_sel[:, 3] == 4)
 
        self.S1T1=np.zeros((56,trial_FR.shape[2]))
        self.S1T2=np.zeros((56,trial_FR.shape[2]))
        self.S2T1=np.zeros((56,trial_FR.shape[2]))
        self.S2T2=np.zeros((56,trial_FR.shape[2]))
        
        
        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(FR_scale[:, :, su_idx]).T

            (base_mean, base_std) = baselineVector(onesu)
            
            onesu=(onesu-base_mean)/base_std
            self.S1T1[:,su_idx] = np.mean(onesu[trial_S1T1,:],axis=0)
            self.S1T2[:,su_idx] = np.mean(onesu[trial_S1T2,:],axis=0)
            self.S2T1[:,su_idx] = np.mean(onesu[trial_S2T1,:],axis=0)
            self.S2T2[:,su_idx] = np.mean(onesu[trial_S2T2,:],axis=0)

            
    def getFeatures(self):
        
        return np.concatenate((self.S1T1,
                               self.S1T2,
                               self.S2T1,
                               self.S2T2))