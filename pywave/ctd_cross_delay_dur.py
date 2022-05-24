# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 09:40:42 2021

@author: Libra
"""


def ctd_cross_delay_dur(feat_per_su, n_neuron=300, n_trial=(20, 25), template_delay=6, bin_range=None, decoder='MCC'):
    template_keys = ["S1_3", "S2_3"] if template_delay == 3 else ["S1_6", "S2_6"]
    score_keys = ["S1_6", "S2_6"] if template_delay == 3 else ["S1_3", "S2_3"]
    avail_sel = [(x[template_keys[0]].shape[1] >= n_trial[1] and x[template_keys[1]].shape[1] >= n_trial[1]
                  and x[score_keys[0]].shape[1] >= n_trial[1] and x[score_keys[1]].shape[1] >= n_trial[1]) for x in
                 feat_per_su]

    if sum(avail_sel) < n_neuron:
        print(f'Not enough SU with suffcient trials {sum(avail_sel)}/{n_neuron}')
        return None

    # bins, trials
    if bin_range is None:
        bin_range = np.arange(feat_per_su[0][template_keys[0]].shape[0])

    scaler = MinMaxScaler()
    if decoder == 'MCC':
        clf = MaximumCorrelationClassifier(n_neuron)
    elif decoder == 'SVC':
        clf = SVC(kernel='linear')

    kf = KFold(10)
    one_cv = []
    for i in range(20):
        su_index = np.random.choice(np.nonzero(avail_sel)[0], n_neuron, replace=False)
        su_selected_features = [feat_per_su[i] for i in su_index]
        template_X1 = []
        score_X1 = []
        template_X2 = []
        score_X2 = []
        template_trial_to_select = np.random.choice(n_trial[1], n_trial[0], replace=False)
        score_trial_to_select = np.random.choice(n_trial[1], n_trial[0], replace=False)
        for one_su in su_selected_features:
            template_X1.append([one_su[template_keys[0]][bin_range, t] for t in template_trial_to_select])
            score_X1.append([one_su[score_keys[0]][bin_range, t] for t in score_trial_to_select])
            template_X2.append([one_su[template_keys[1]][bin_range, t] for t in template_trial_to_select])
            score_X2.append([one_su[score_keys[1]][bin_range, t] for t in score_trial_to_select])

        template_X1 = np.array(template_X1).transpose((1, 0, 2))
        template_X2 = np.array(template_X2).transpose((1, 0, 2))  # trial, SU, bin
        score_X1 = np.array(score_X1).transpose((1, 0, 2))
        score_X2 = np.array(score_X2).transpose((1, 0, 2))  # trial, SU, bin

        for (templates, tests) in kf.split(template_X1):
            score_mat = np.zeros((bin_range.shape[0], bin_range.shape[0]))
            for template_bin_idx in np.arange(bin_range.shape[0]):
                X_templates = np.vstack(
                    (template_X1[templates, :, template_bin_idx], template_X2[templates, :, template_bin_idx]))
                y_templates = np.hstack((np.zeros_like(templates), np.ones_like(templates))).T
                scaler = scaler.fit(X_templates)
                X_templates = scaler.transform(X_templates)
                clf.fit(X_templates, y_templates)
                for test_bin_idx in np.arange(bin_range.shape[0]):
                    X_test = np.vstack((score_X1[tests, :, test_bin_idx], score_X2[tests, :, test_bin_idx]))
                    y_test = np.hstack((np.zeros_like(tests), np.ones_like(tests))).T
                    # y_shuf=y.copy()
                    # rng.shuffle(y_shuf)
                    X_test = scaler.transform(X_test)
                    if decoder == 'MCC':
                        score_mat[template_bin_idx, test_bin_idx] = np.mean(np.hstack(
                            (clf.predict(score_X1[tests, :, test_bin_idx]) <= 0,
                             clf.predict(score_X2[tests, :, test_bin_idx]) > 0)))
                    elif decoder == 'SVC':
                        score_mat[template_bin_idx, test_bin_idx] = clf.score(X_test, y_test)

            score_mat = score_mat * 100
            one_cv.append(score_mat)

    return one_cv

    ### disabled due to missing trials
    # availErrTrials=[]
    # for su in features_per_su:
    #     availErrTrials.append([su['S1_3_ERR'].shape[1],su['S2_3_ERR'].shape[1],su['S1_6_ERR'].shape[1],su['S2_6_ERR'].shape[1]])