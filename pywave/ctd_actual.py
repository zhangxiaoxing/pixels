# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 09:44:33 2021

@author: Libra
"""

import numpy as np

from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVC
from mcc.mcc import MaximumCorrelationClassifier

def ctd_actual(fr_per_su, n_neuron=50, n_trial=(20, 25), delay=6, bin_range=None, decoder='MCC'):
    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]
    avail_sel = [(x[keys[0]].shape[1] >= n_trial[1] and x[keys[1]].shape[1] >= n_trial[1]) for x in fr_per_su]

    if sum(avail_sel) < n_neuron:
        print('Not enough SU with suffcient trials')
        return None

    # bins, trials
    if bin_range is None:
        bin_range = np.arange(fr_per_su[0][keys[0]].shape[0])

    scaler = MinMaxScaler()
    if decoder == 'MCC':
        clf = MaximumCorrelationClassifier(n_neuron)
    elif decoder == 'SVC':
        clf = SVC(kernel='linear')

    kf = KFold(10)
    one_cv = []
    for i in range(20):
        su_index = np.random.choice(np.nonzero(avail_sel)[0], n_neuron, replace=False)
        su_selected_features = [fr_per_su[i] for i in su_index]
        X1 = []
        X2 = []
        trial_to_select = np.random.choice(n_trial[1], n_trial[0], replace=False)
        for one_su in su_selected_features:
            X1.append([one_su[keys[0]][bin_range, t] for t in trial_to_select])
            X2.append([one_su[keys[1]][bin_range, t] for t in trial_to_select])

        X1 = np.array(X1).transpose((1, 0, 2))
        X2 = np.array(X2).transpose((1, 0, 2))  # trial, SU, bin

        for (templates, tests) in kf.split(X1):
            score_mat = np.ones((bin_range.shape[0], bin_range.shape[0]))
            for template_bin_idx in np.arange(bin_range.shape[0]):
                X_templates = np.vstack((X1[templates, :, template_bin_idx], X2[templates, :, template_bin_idx]))
                y_templates = np.hstack((np.ones_like(templates) * -1, np.ones_like(templates) * 1)).T
                # np.random.shuffle(y_templates)
                scaler = scaler.fit(X_templates)
                X_templates = scaler.transform(X_templates)
                clf.fit(X_templates, y_templates)
                for test_bin_idx in np.arange(bin_range.shape[0]):
                    X_test = np.vstack((X1[tests, :, test_bin_idx], X2[tests, :, test_bin_idx]))
                    y_test = np.hstack((np.ones_like(tests) * -1, np.ones_like(tests) *1)).T
                    # y_shuf=y.copy()
                    # rng.shuffle(y_shuf)
                    X_test = scaler.transform(X_test)
                    if decoder == 'MCC':
                        score_mat[template_bin_idx, test_bin_idx] = np.mean(np.hstack(
                            (clf.predict(X1[tests, :, test_bin_idx]) <= 0,
                             clf.predict(X2[tests, :, test_bin_idx]) > 0)))
                    else:
                        score_mat[template_bin_idx, test_bin_idx] = clf.score(X_test, y_test)

            score_mat = score_mat * 100
            one_cv.append(score_mat)

    return one_cv