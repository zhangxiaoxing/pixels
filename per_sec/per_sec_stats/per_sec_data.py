# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 11:31:50 2021

@author: Libra
"""
import numpy as np

class per_sec_data:
    def __init__(self):
        self.per_sec_fr = None # 28 x su, S1, S2
        self.per_sec_sel = None  # 14 x SU #selectivity
        self.per_sec_auc = None  # 14 x SU
        self.per_sec_wrs_p = None  # 7 x SU , sample + 6s delay

    def getFeatures(self):
        return (self.per_sec_sel, self.non_sel_mod, self.per_sec_prefS1, self.per_sec_prefS2, self.baseline_sel,
                self.per_sec_sel_raw,self.per_sec_auc,self.per_sec_wrs_p,self.per_sec_fr)