# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 10:28:51 2021

@author: Libra
"""

import numpy as np
import scipy.stats as stats
from sklearn.metrics import roc_auc_score


def wrs(A, B, bonf=1):
####TODO: fischers exact alternative
####TODO: check for variance


    (stat, p) = stats.mannwhitneyu(
        A.flatten(),
        B.flatten(),
        alternative="two-sided")
    return p*bonf

def auc(A,B):
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

    return auc


def exact_mc_perm_test(self, xs, ys, nmc):
    n, k = len(xs), 0
    diff = np.abs(np.mean(xs) - np.mean(ys))
    zs = np.concatenate([xs, ys])
    for j in range(nmc):
        np.random.shuffle(zs)
        k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
    return k / nmc