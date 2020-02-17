# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 08:20:49 2020

@author: Libra
"""

import glob
import numpy as np
import matplotlib.pyplot as plt


params=np.load('consec_nonneg_trials160_opti_params.npy')
plt.figure()
plt.plot(params[:,0],params[:,1])
plt.plot(params[:,0],params[:,3])


fl=glob.glob('*.npz')
for fpath in fl:
    fstr=np.load(fpath)
    