# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 10:28:51 2021

@author: Libra
"""

import numpy as np
import scipy.stats as stats
from sklearn.metrics import roc_auc_score


def bool_stats_test(self, A, B, bonf=1):
####TODO: fischers exact alternative
####TODO: check for variance

    try:
        # flatten is used instead of mean to control for variance
        # ranksum is supposedly insensitive to length of data
        (stat, p) = stats.mannwhitneyu(
            A.flatten(),
            B.flatten(),
            alternative="two-sided",
        )
    except ValueError:
        p = 1

    try:
        auc = roc_auc_score(
            np.concatenate((
                np.zeros(A.size),
                np.ones(B.size)
            )),
            np.concatenate((
                A.flatten(),
                B.flatten(),
            ))
        )
    except ValueError:
        auc = 0.5

    if np.sum(A)+np.sum(B)==0:
        selectivity=0
    else:
        mma = np.mean(A)
        mmb = np.mean(B)
        selectivity=(mma-mmb)/(mma+mmb)

    return (p < (0.05 / bonf), p, selectivity,auc)


def exact_mc_perm_test(self, xs, ys, nmc):
    n, k = len(xs), 0
    diff = np.abs(np.mean(xs) - np.mean(ys))
    zs = np.concatenate([xs, ys])
    for j in range(nmc):
        np.random.shuffle(zs)
        k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
    return k / nmc