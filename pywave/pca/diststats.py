# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:22:35 2021

@author: Libra
"""

import numpy as np
from sklearn.decomposition import PCA

from scipy.spatial.distance import euclidean


class Dist_stats:
    def __init__(self, ref_comp):
        self.ref_comp = ref_comp
        self.S1_S2_3_dist = []
        self.S1_S2_6_dist = []
        self.exit_dist_3_S1 = []
        self.exit_dist_3_S2 = []
        self.exit_dist_6_S1 = []
        self.exit_dist_6_S2 = []

        self.latedelay_dist_3_S1 = []
        self.latedelay_dist_3_S2 = []
        self.latedelay_dist_6_S1 = []
        self.latedelay_dist_6_S2 = []

        self.ref_dist_3_S1 = []
        self.ref_dist_3_S2 = []
        self.ref_dist_6_S1 = []
        self.ref_dist_6_S2 = []

        self.coeff = []

    def append_data(self, fr, coeff):
        pcamat = np.hstack(
            [
                np.vstack([np.mean(
                    np.vstack((x["S1_S1_3m"][8:28], x["S1_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S2_S1_3m"][8:28], x["S2_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S1_S1_6m"][8:40], x["S1_S2_6m"][8:40])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S2_S1_6m"][8:40], x["S2_S2_6m"][8:40])), axis=0) for x in fr]),
            ]
        )
        # pca = PCA(n_components=20)
        # comp = pca.fit_transform(pcamat.T)
        # coeff = pca.components_
        # ratio = pca.explained_variance_ratio_

        pcamat_cent = pcamat - np.expand_dims(np.mean(pcamat, axis=1), axis=1)
        comp = np.matmul(coeff, pcamat_cent).T

        S1_S2_3 = []
        deci_vec_S1_3 = []
        deci_vec_S2_3 = []
        late_vec_S1_3 = []
        late_vec_S2_3 = []
        ref_S1_3 = []
        ref_S2_3 = []
        for bin in range(20):  # time range
            # comp range
            S1_S2_3.append(euclidean(comp[bin, :3], comp[bin + 20, :3]))
            deci_vec_S1_3.append(
                euclidean(comp[bin, :3], self.ref_comp[19, :3]))
            deci_vec_S2_3.append(
                euclidean(comp[bin + 20, :3], self.ref_comp[39, :3]))
            late_vec_S1_3.append(
                euclidean(comp[bin, :3], self.ref_comp[17, :3]))
            late_vec_S2_3.append(
                euclidean(comp[bin + 20, :3], self.ref_comp[37, :3]))
            ref_S1_3.append(euclidean(comp[bin, :3], self.ref_comp[bin, :3]))
            ref_S2_3.append(
                euclidean(comp[bin + 20, :3], self.ref_comp[bin + 20, :3]))
        self.S1_S2_3_dist.append(S1_S2_3)
        self.exit_dist_3_S1.append(deci_vec_S1_3)
        self.exit_dist_3_S2.append(deci_vec_S2_3)
        self.latedelay_dist_3_S1.append(late_vec_S1_3)
        self.latedelay_dist_3_S2.append(late_vec_S2_3)
        self.ref_dist_3_S1.append(ref_S1_3)
        self.ref_dist_3_S2.append(ref_S2_3)
        S1_S2_6 = []
        deci_vec_S1_6 = []
        deci_vec_S2_6 = []
        late_vec_S1_6 = []
        late_vec_S2_6 = []
        ref_S1_6 = []
        ref_S2_6 = []
        for bin in range(40, 72):  # time range
            # comp range
            S1_S2_6.append(euclidean(comp[bin, :3], comp[bin + 32, :3]))
            deci_vec_S1_6.append(
                euclidean(comp[bin, :3], self.ref_comp[19, :3]))
            deci_vec_S2_6.append(
                euclidean(comp[bin + 32, :3], self.ref_comp[39, :3]))
            late_vec_S1_6.append(
                euclidean(comp[bin, :3], self.ref_comp[17, :3]))
            late_vec_S2_6.append(
                euclidean(comp[bin + 32, :3], self.ref_comp[37, :3]))
            ref_S1_6.append(euclidean(comp[bin, :3], self.ref_comp[bin, :3]))
            ref_S2_6.append(
                euclidean(comp[bin + 32, :3], self.ref_comp[bin + 32, :3]))

        self.S1_S2_6_dist.append(S1_S2_6)
        self.exit_dist_6_S1.append(deci_vec_S1_6)
        self.exit_dist_6_S2.append(deci_vec_S2_6)
        self.latedelay_dist_6_S1.append(late_vec_S1_6)
        self.latedelay_dist_6_S2.append(late_vec_S2_6)
        self.ref_dist_6_S1.append(ref_S1_6)
        self.ref_dist_6_S2.append(ref_S2_6)

    def get_coeff(self, fr):
        pcamat = np.hstack(
            [
                np.vstack([np.mean(
                    np.vstack((x["S1_S1_3m"][8:28], x["S1_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S2_S1_3m"][8:28], x["S2_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S1_S1_6m"][8:40], x["S1_S2_6m"][8:40])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S2_S1_6m"][8:40], x["S2_S2_6m"][8:40])), axis=0) for x in fr]),
            ]
        )
        pca = PCA(n_components=20)
        comp = pca.fit_transform(pcamat.T)
        coeff = pca.components_
        # ratio = pca.explained_variance_ratio_

        return coeff

    def get_data(self):
        return (self.S1_S2_3_dist, self.S1_S2_6_dist, self.exit_dist_3_S1,
                self.exit_dist_3_S2, self.exit_dist_6_S1, self.exit_dist_6_S2,
                self.ref_dist_3_S1, self.ref_dist_3_S2, self.ref_dist_6_S1, self.ref_dist_6_S2,
                self.latedelay_dist_3_S1, self.latedelay_dist_3_S2, self.latedelay_dist_6_S1, self.latedelay_dist_6_S2)