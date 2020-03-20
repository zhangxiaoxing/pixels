"""
This code calculates cross time decoding and cross delay duration decoding, generates heatmaps
depends on the per_second_stats module for sustained and transient classification


"""

import sys
import per_second_stats
import os
import numpy as np

from sklearn.svm import LinearSVC
from sklearn.svm import SVC
# from mcc.mcc import MaximumCorrelationClassifier


from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler
from multiprocessing import Pool
import matplotlib.pyplot as plt
import scipy.stats as stats


def exact_mc_perm_test(xs, ys, nmc):
    n, k = len(xs), 0
    diff = np.abs(np.mean(xs) - np.mean(ys))
    zs = np.concatenate([xs, ys])
    for j in range(nmc):
        np.random.shuffle(zs)
        k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
    return k / nmc


def last_bin_decoding(delay=6):
    fstr = np.load("ctd.npz", allow_pickle=True)
    features_per_su = fstr["features_per_su"].tolist()
    baseline_WRS_p3_p6 = baseline_statstics(features_per_su)
    wrs_p = baseline_WRS_p3_p6[:, 0] if delay == 3 else baseline_WRS_p3_p6[:, 1]
    # reg_list = fstr["reg_list"].tolist()
    sus_trans_flag = per_second_stats.process_all()  # 33172 x 4, sust,trans,switch,unclassified
    sus_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[0, :], wrs_p > 0.05))[0]]
    trans_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[1, :], wrs_p > 0.05))[0]]

    sust_sum = []
    repeats = 1000
    for n_neuron in np.arange(25, 126, 25):
        sust_sum.append(same_time_decoding(sus_feat, n_neuron, (20, 25), 6, np.arange(36, 40), repeats))
    trans_sum = []
    for n_neuron in np.arange(25, 2001, 25):
        trans_sum.append(same_time_decoding(trans_feat, n_neuron, (20, 25), 6, np.arange(36, 40), repeats))
    return (sust_sum, trans_sum)


def same_time_decoding(features_per_su, n_neuron=300, n_trial=(20, 25), delay=6, bin_range=None, repeat=10):
    if bin_range is None:
        bin_range = np.arange(36, 40)

    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]
    concurrent_time = []
    avail_sel = [(x[keys[0]].shape[1] >= n_trial[1] and x[keys[1]].shape[1] >= n_trial[1]) for x in features_per_su]
    for i in range(repeat):
        su_index = np.random.choice(np.nonzero(avail_sel)[0], n_neuron, replace=False)
        su_selected_features = [features_per_su[i] for i in su_index]
        X1 = []
        X2 = []
        trial_to_select = np.random.choice(n_trial[1], n_trial[0], replace=False)
        for one_su in su_selected_features:
            X1.append([np.mean(one_su[keys[0]][bin_range, t]) for t in trial_to_select])
            X2.append([np.mean(one_su[keys[1]][bin_range, t]) for t in trial_to_select])

        # bins, trials
        rng = np.random.default_rng()
        X1 = np.array(X1).T
        X2 = np.array(X2).T
        y1 = np.ones((X1.shape[0]))
        y2 = np.zeros_like(y1)
        y = np.hstack((y1, y2))
        y_shuf = y.copy()
        rng.shuffle(y_shuf)
        scaler = MinMaxScaler()
        X = scaler.fit_transform(np.vstack((X1, X2)))
        clf = LinearSVC()
        scores = cross_val_score(clf, X, y, cv=10, n_jobs=-1) * 100
        scores_shuffled = cross_val_score(clf, X, y_shuf, cv=10, n_jobs=-1) * 100
        concurrent_time.append(
            [
                np.mean(scores),
                np.std(scores),
                np.mean(scores_shuffled),
                np.std(scores_shuffled),
            ]
        )
    return concurrent_time


def cross_time_decoding(denovo=False, to_plot=False, delay=6, cpu=30):
    repeats = 1000
    sus50 = None
    trans50 = None
    trans1000 = None
    if denovo:
        fstr = np.load("ctd.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        baseline_WRS_p3_p6 = baseline_statstics(features_per_su)
        wrs_p = baseline_WRS_p3_p6[:, 0] if delay == 3 else baseline_WRS_p3_p6[:, 1]
        # reg_list = fstr["reg_list"].tolist()
        sus_trans_flag = per_second_stats.process_all()  # 33172 x 4, sust,trans,switch,unclassified
        sus_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[0, :], wrs_p > 0.05))[0]]
        trans_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[1, :], wrs_p > 0.05))[0]]
        curr_pool = Pool(processes=cpu)
        sus50 = []
        trans50 = []
        trans1000 = []
        sust_proc = []
        trans50_proc = []
        trans1000_proc = []
        for i in range(repeats):
            sust_proc.append(curr_pool.apply_async(cross_time_decoding_actual, args=(sus_feat,),
                                                   kwds={"n_neuron": 50, "n_trial": (20, 25), "delay": delay,
                                                         "bin_range": np.arange(4, 52)}))

            trans50_proc.append(curr_pool.apply_async(cross_time_decoding_actual, args=(trans_feat,),
                                                      kwds={"n_neuron": 50, "n_trial": (20, 25), "delay": delay,
                                                            "bin_range": np.arange(4, 52)}))
            trans1000_proc.append(curr_pool.apply_async(cross_time_decoding_actual, args=(trans_feat,),
                                                        kwds={"n_neuron": 1000, "n_trial": (20, 25), "delay": delay,
                                                              "bin_range": np.arange(4, 52)}))
        for one_proc in sust_proc:
            sus50.append(one_proc.get())

        for one_proc in trans50_proc:
            trans50.append(one_proc.get())

        for one_proc in trans1000_proc:
            trans1000.append(one_proc.get())

        curr_pool.close()
        curr_pool.join()
        np.savez_compressed(f'sus_trans_ctd_{delay}_{repeats}.npz', sus50=sus50, trans50=trans50,
                            trans1000=trans1000)
    else:
        fstr = np.load(os.path.join('ctd', f'sus_trans_ctd_{delay}_1000.npz'), 'r')
        sus50 = fstr['sus50']
        trans50 = fstr['trans50']
        trans1000 = fstr['trans1000']

    if to_plot:
        (fig, ax) = plt.subplots(1, 3, figsize=[9, 3], dpi=200)
        ax[0].imshow(sus50.mean(axis=0).mean(axis=0), cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[0].set_title('50 sustained SU')
        ax[0].set_ylabel('template time (s)')
        ax[1].imshow(trans50.mean(axis=0).mean(axis=0), cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[1].set_title('50 transient SU')
        im2 = ax[2].imshow(trans1000.mean(axis=0).mean(axis=0), cmap="jet", aspect="auto", origin='lower', vmin=0,
                           vmax=100)
        ax[2].set_title('1000 transient SU')
        plt.colorbar(im2, ticks=[0, 50, 100], format="%d")

        for oneax in ax:
            oneax.set_xticks([7.5, 27.5])
            oneax.set_xticklabels([0, 5])
            oneax.set_xlabel('scoring time (s)')
            oneax.set_yticks([7.5, 27.5])
            oneax.set_yticklabels([0, 5])

            if delay == 6:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
            elif delay == 3:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]

        fig.savefig(f'ctd_{delay}_{repeats}.png', bbox_inches='tight')

        plt.show()
        return (sus50, trans50, trans1000)


def ctd_cross_all(denovo=False, to_plot=False, delay=3, cpu=30):
    repeats = 1000
    if denovo:
        fstr = np.load("ctd.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        baseline_WRS_p3_p6 = baseline_statstics(features_per_su)
        wrs_p = baseline_WRS_p3_p6[:, 0] if delay == 3 else baseline_WRS_p3_p6[:, 1]
        # reg_list = fstr["reg_list"].tolist()
        sus_trans_flag = per_second_stats.process_all()  # 33172 x 4, sust,trans,switch,unclassified
        sus_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[0, :], wrs_p > 0.05))[0]]
        trans_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[1, :], wrs_p > 0.05))[0]]
        curr_pool = Pool(processes=cpu)
        sus50 = []
        trans50 = []
        trans1000 = []
        sust_proc = []
        trans50_proc = []
        trans1000_proc = []
        for i in range(repeats):
            sust_proc.append(curr_pool.apply_async(cross_time_decoding_cross, args=(sus_feat,),
                                                   kwds={"n_neuron": 50, "n_trial": (20, 25), "template_delay": delay,
                                                         "bin_range": np.arange(4, 52)}))

            trans50_proc.append(curr_pool.apply_async(cross_time_decoding_cross, args=(trans_feat,),
                                                      kwds={"n_neuron": 50, "n_trial": (20, 25),
                                                            "template_delay": delay,
                                                            "bin_range": np.arange(4, 52)}))
            trans1000_proc.append(curr_pool.apply_async(cross_time_decoding_cross, args=(trans_feat,),
                                                        kwds={"n_neuron": 1000, "n_trial": (20, 25),
                                                              "template_delay": delay,
                                                              "bin_range": np.arange(4, 52)}))
        for one_proc in sust_proc:
            sus50.append(one_proc.get())

        for one_proc in trans50_proc:
            trans50.append(one_proc.get())

        for one_proc in trans1000_proc:
            trans1000.append(one_proc.get())

        curr_pool.close()
        curr_pool.join()
        np.savez_compressed(f'sus_trans_ctd_cross_delay_dur_{delay}_{repeats}.npz', sus50=sus50, trans50=trans50,
                            trans1000=trans1000)
    else:
        fstr = np.load(f"sus_trans_ctd_{delay}_{repeats}.npz")
        sus50 = fstr['sus50']
        trans50 = fstr['trans50']
        trans1000 = fstr['trans1000']

    sus_mat = np.array(sus50).mean(axis=0).mean(axis=0)
    trans50_mat = np.array(trans50).mean(axis=0).mean(axis=0)
    trans1000_mat = np.array(trans1000).mean(axis=0).mean(axis=0)
    if to_plot:
        (fig, ax) = plt.subplots(1, 3, figsize=[9, 3], dpi=200)
        ax[0].imshow(sus_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[0].set_title('50 sustained SU')
        ax[0].set_ylabel('template time (s)')
        ax[1].imshow(trans50_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[1].set_title('50 transient SU')
        im2 = ax[2].imshow(trans1000_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[2].set_title('1000 transient SU')
        plt.colorbar(im2, ticks=[0, 50, 100], format="%d")

        for oneax in ax:
            oneax.set_xticks([7.5, 27.5])
            oneax.set_xticklabels([0, 5])
            oneax.set_xlabel('scoring time (s)')
            oneax.set_yticks([7.5, 27.5])
            oneax.set_yticklabels([0, 5])

            if delay == 6:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]
            elif delay == 3:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]

        fig.savefig('ctd_cross_delay_dur.png', bbox_inches='tight')
        plt.show()

    return (sus50, trans50, trans1000)


def ctd_correct_error_all(denovo=False, to_plot=False, delay=3, cpu=30):
    repeats = 1000
    if denovo:
        fstr = np.load("ctd.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        baseline_WRS_p3_p6 = baseline_statstics(features_per_su)
        wrs_p = baseline_WRS_p3_p6[:, 0] if delay == 3 else baseline_WRS_p3_p6[:, 1]
        # reg_list = fstr["reg_list"].tolist()
        sus_trans_flag = per_second_stats.process_all()  # 33172 x 4, sust,trans,switch,unclassified
        sus_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[0, :], wrs_p > 0.05))[0]]
        trans_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[1, :], wrs_p > 0.05))[0]]
        curr_pool = Pool(processes=cpu)
        sus50 = []
        trans50 = []
        trans1000 = []
        sust_proc = []
        trans50_proc = []
        trans1000_proc = []
        for i in range(repeats):
            sust_proc.append(curr_pool.apply_async(ctd_correct_error, args=(sus_feat,),
                                                   kwds={"n_neuron": 50, "n_trial": (20, 25, 2, 4),
                                                         "template_delay": delay,
                                                         "bin_range": np.arange(4, 52)}))

            trans50_proc.append(curr_pool.apply_async(ctd_correct_error, args=(trans_feat,),
                                                      kwds={"n_neuron": 50, "n_trial": (20, 25, 2, 4),
                                                            "template_delay": delay,
                                                            "bin_range": np.arange(4, 52)}))
            trans1000_proc.append(curr_pool.apply_async(ctd_correct_error, args=(trans_feat,),
                                                        kwds={"n_neuron": 1000, "n_trial": (20, 25, 2, 4),
                                                              "template_delay": delay,
                                                              "bin_range": np.arange(4, 52)}))
        for one_proc in sust_proc:
            sus50.append(one_proc.get())

        for one_proc in trans50_proc:
            trans50.append(one_proc.get())

        for one_proc in trans1000_proc:
            trans1000.append(one_proc.get())

        curr_pool.close()
        curr_pool.join()
        np.savez_compressed(f'sus_trans_ctd_correct_error_{delay}_{repeats}.npz', sus50=sus50, trans50=trans50,
                            trans1000=trans1000)
    else:
        fstr = np.load(f'sus_trans_ctd_correct_error_{delay}_{repeats}.npz')
        sus50 = fstr['sus50']
        trans50 = fstr['trans100']
        trans1000 = fstr['trans1000']

    sus_mat = np.array(sus50).mean(axis=0).mean(axis=0)
    trans50_mat = np.array(trans50).mean(axis=0).mean(axis=0)
    trans1000_mat = np.array(trans1000).mean(axis=0).mean(axis=0)
    if to_plot:
        (fig, ax) = plt.subplots(1, 3, figsize=[9, 3], dpi=200)
        ax[0].imshow(sus_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[0].set_title('50 sustained SU')
        ax[0].set_ylabel('template time (s)')
        ax[1].imshow(trans50_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[1].set_title('50 transient SU')
        im2 = ax[2].imshow(trans1000_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[2].set_title('1000 transient SU')
        plt.colorbar(im2, ticks=[0, 50, 100], format="%d")

        for oneax in ax:
            oneax.set_xticks([7.5, 27.5])
            oneax.set_xticklabels([0, 5])
            oneax.set_xlabel('scoring time (s)')
            oneax.set_yticks([7.5, 27.5])
            oneax.set_yticklabels([0, 5])

            if delay == 6:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
            elif delay == 3:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]

        fig.savefig(f'ctd_correct_error_{delay}.png', bbox_inches='tight')
        plt.show()

    return (sus50, trans50, trans1000)


def cross_time_decoding_actual(features_per_su, n_neuron=50, n_trial=(20, 25), delay=6, bin_range=None):
    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]
    avail_sel = [(x[keys[0]].shape[1] >= n_trial[1] and x[keys[1]].shape[1] >= n_trial[1]) for x in features_per_su]

    if sum(avail_sel) < n_neuron:
        print('Not enough SU with suffcient trials')
        return None

    # bins, trials
    if bin_range is None:
        bin_range = np.arange(features_per_su[0][keys[0]].shape[0])

    scaler = MinMaxScaler()
    clf = SVC(kernel='linear')
    # clf=MaximumCorrelationClassifier(n_neuron)
    kf = KFold(10)
    su_index = np.random.choice(np.nonzero(avail_sel)[0], n_neuron, replace=False)
    su_selected_features = [features_per_su[i] for i in su_index]
    X1 = []
    X2 = []
    trial_to_select = np.random.choice(n_trial[1], n_trial[0], replace=False)
    for one_su in su_selected_features:
        X1.append([one_su[keys[0]][bin_range, t] for t in trial_to_select])
        X2.append([one_su[keys[1]][bin_range, t] for t in trial_to_select])

    X1 = np.array(X1).transpose((1, 0, 2))
    X2 = np.array(X2).transpose((1, 0, 2))  # trial, SU, bin

    one_cv = []
    for (templates, tests) in kf.split(X1):
        score_mat = np.ones((bin_range.shape[0], bin_range.shape[0]))
        for template_bin_idx in np.arange(bin_range.shape[0]):
            X_templates = np.vstack((X1[templates, :, template_bin_idx], X2[templates, :, template_bin_idx]))
            y_templates = np.hstack((np.ones_like(templates) * 1, np.ones_like(templates) * 2)).T
            # np.random.shuffle(y_templates)
            scaler = scaler.fit(X_templates)
            X_templates = scaler.transform(X_templates)
            for test_bin_idx in np.arange(bin_range.shape[0]):
                X_test = np.vstack((X1[tests, :, test_bin_idx], X2[tests, :, test_bin_idx]))
                y_test = np.hstack((np.ones_like(tests) * 1, np.ones_like(tests) * 2)).T
                # y_shuf=y.copy()
                # rng.shuffle(y_shuf)
                clf.fit(X_templates, y_templates)
                X_test = scaler.transform(X_test)
                # score_mat[template_bin_idx, test_bin_idx] = np.mean(np.hstack((clf.predict(X1[tests, :, test_bin_idx])<=0,clf.predict(X2[tests, :, test_bin_idx])>0)))

                score_mat[template_bin_idx, test_bin_idx] = clf.score(X_test, y_test)

        score_mat = score_mat * 100
        one_cv.append(score_mat)
    np.mean(one_cv)

    return one_cv


def cross_time_decoding_cross(features_per_su, n_neuron=300, n_trial=(20, 25), template_delay=6, bin_range=None):
    template_keys = ["S1_3", "S2_3"] if template_delay == 3 else ["S1_6", "S2_6"]
    score_keys = ["S1_6", "S2_6"] if template_delay == 3 else ["S1_3", "S2_3"]
    avail_sel = [(x[template_keys[0]].shape[1] >= n_trial[1] and x[template_keys[1]].shape[1] >= n_trial[1]
                  and x[score_keys[0]].shape[1] >= n_trial[1] and x[score_keys[1]].shape[1] >= n_trial[1]) for x in
                 features_per_su]

    if sum(avail_sel) < n_neuron:
        print('Not enough SU with suffcient trials')
        return None

    # bins, trials
    if bin_range is None:
        bin_range = np.arange(features_per_su[0][template_keys[0]].shape[0])

    scaler = MinMaxScaler()
    clf = LinearSVC()
    kf = KFold(10)
    su_index = np.random.choice(np.nonzero(avail_sel)[0], n_neuron, replace=False)
    su_selected_features = [features_per_su[i] for i in su_index]
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

    one_cv = []
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
                score_mat[template_bin_idx, test_bin_idx] = clf.score(X_test, y_test)

        score_mat = score_mat * 100
        one_cv.append(score_mat)

    return one_cv

    ### disabled due to missing trials
    # availErrTrials=[]
    # for su in features_per_su:
    #     availErrTrials.append([su['S1_3_ERR'].shape[1],su['S2_3_ERR'].shape[1],su['S1_6_ERR'].shape[1],su['S2_6_ERR'].shape[1]])


def ctd_correct_error(features_per_su, n_neuron=300, n_trial=(20, 25, 2, 4), template_delay=6, bin_range=None):
    template_keys = ["S1_3", "S2_3"] if template_delay == 3 else ["S1_6", "S2_6"]
    score_keys = ["S1_3_ERR", "S2_3_ERR"] if template_delay == 3 else ["S1_6_ERR", "S2_6_ERR"]
    avail_sel = [((x[score_keys[0]].ndim > 1) and (x[score_keys[0]].shape[1] > n_trial[3]) and (
            x[score_keys[1]].ndim > 1) and (x[score_keys[1]].shape[1] > n_trial[3]) and (
                          x[template_keys[0]].shape[1] >= n_trial[1]) and (x[template_keys[1]].shape[1] >= n_trial[1]
                                                                           )) for x in features_per_su]

    if sum(avail_sel) < n_neuron:
        print('Not enough SU with suffcient trials')
        return None

    # bins, trials
    if bin_range is None:
        bin_range = np.arange(features_per_su[0][template_keys[0]].shape[0])

    scaler = MinMaxScaler()
    clf = LinearSVC()
    kf = KFold(10)
    su_index = np.random.choice(np.nonzero(avail_sel)[0], n_neuron, replace=False)
    su_selected_features = [features_per_su[i] for i in su_index]
    template_X1 = []
    score_X1 = []
    template_X2 = []
    score_X2 = []
    template_trial_to_select = np.random.choice(n_trial[1], n_trial[0], replace=False)
    score_trial_to_select = np.random.choice(n_trial[3], n_trial[2], replace=False)
    for one_su in su_selected_features:
        template_X1.append([one_su[template_keys[0]][bin_range, t] for t in template_trial_to_select])
        score_X1.append([one_su[score_keys[0]][bin_range, t] for t in score_trial_to_select])
        template_X2.append([one_su[template_keys[1]][bin_range, t] for t in template_trial_to_select])
        score_X2.append([one_su[score_keys[1]][bin_range, t] for t in score_trial_to_select])

    template_X1 = np.array(template_X1).transpose((1, 0, 2))
    template_X2 = np.array(template_X2).transpose((1, 0, 2))  # trial, SU, bin
    score_X1 = np.array(score_X1).transpose((1, 0, 2))
    score_X2 = np.array(score_X2).transpose((1, 0, 2))  # trial, SU, bin

    one_cv = []
    (templates, tests) = list(kf.split(template_X1))[0]

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
            score_mat[template_bin_idx, test_bin_idx] = clf.score(X_test, y_test)

    score_mat = score_mat * 100
    one_cv.append(score_mat)

    return one_cv

    ### disabled due to missing trials
    # availErrTrials=[]
    # for su in features_per_su:
    #     availErrTrials.append([su['S1_3_ERR'].shape[1],su['S2_3_ERR'].shape[1],su['S1_6_ERR'].shape[1],su['S2_6_ERR'].shape[1]])


def statistical_test_correct_error(delay=3):
    fstr = np.load(os.path.join('ctd', f'sus_trans_ctd_{delay}_1000.npz'), 'r')
    fstr_Err = np.load(os.path.join('ctd', f'sus_trans_ctd_correct_error_{delay}_1000.npz'), 'r')

    correct_mat = fstr['trans1000']
    err_mat = fstr_Err['trans1000']

    err_mean = np.ones((correct_mat.shape[2:4]))
    cr_mean = np.ones((correct_mat.shape[2:4]))
    for template_bin in range(correct_mat.shape[2]):
        for score_bin in range(correct_mat.shape[3]):
            err_mean[template_bin, score_bin] = np.mean(err_mat[:, :, template_bin, score_bin])
            cr_mean[template_bin, score_bin] = np.mean(correct_mat[:, :, template_bin, score_bin])

    p_mat = np.ones((correct_mat.shape[2:4]))
    for template_bin in range(correct_mat.shape[2]):
        for score_bin in range(correct_mat.shape[3]):
            (_stat, p_mat[template_bin, score_bin]) = stats.mannwhitneyu(
                correct_mat[:, :, template_bin, score_bin].flatten(),
                err_mat[:, :, template_bin, score_bin].flatten(),
                alternative="two-sided")
    p_mat = p_mat * correct_mat.shape[2] * correct_mat.shape[3]
    (fig, ax) = plt.subplots(1, 1)
    ax.imshow(p_mat > 0.05)
    plt.show()


def baseline_statstics(features_per_su):
    baseline_stats = []
    for su in features_per_su:
        S1_3 = su['S1_3'][2:10, :].flatten()
        S2_3 = su['S2_3'][2:10, :].flatten()
        S1_6 = su['S1_6'][2:10, :].flatten()
        S2_6 = su['S2_6'][2:10, :].flatten()
        p3 = None
        p6 = None
        if not (np.unique(np.concatenate((S1_3, S2_3))).shape[0] == 1):
            p3 = stats.mannwhitneyu(S1_3, S2_3, alternative='two-sided')[1]
        else:
            p3 = 1

        if not (np.unique(np.concatenate((S1_6, S2_6))).shape[0] == 1):
            p6 = stats.mannwhitneyu(S1_6, S2_6, alternative='two-sided')[1]
        else:
            p6 = 1
        baseline_stats.append((p3, p6))

    np.savez_compressed('baseline_WRS.npz', baseline_WRS_p3_p6=baseline_stats)
    return np.array(baseline_stats)

    # idx=np.nonzero([stat.pvalue < 0.05 for stat in baseline_stats])[0]

    # idx = np.nonzero(np.array(baseline_sel) > 0.99)[0]

    # s1m=np.mean(features_per_su[idx[2]]['S1_6'],axis=1)
    # s2m=np.mean(features_per_su[idx[2]]['S2_6'],axis=1)
    # for i in range(len(idx)):
    #     (fig,ax)=plt.subplots(2,1,sharex=True,figsize=(4,4),dpi=200)
    #     ax[0].imshow(features_per_su[idx[i]]['S1_6'].T,cmap='jet', aspect='auto')
    #     ax[0].set_ylabel('S1 trial #')
    #     [ax[0].axvline(x,color='w',ls=':') for x in [11.5,15.5, 39.5,43.5]]
    #     [ax[0].axvline(x, color='gray', ls=':') for x in [1.5, 9.5]]
    #     ax[1].imshow(features_per_su[idx[i]]['S2_6'].T,cmap='jet', aspect='auto')
    #     ax[1].set_ylabel('S2 trial #')
    #     ax[1].set_xticks([11.5,31.5,51.5])
    #     [ax[1].axvline(x, color='w', ls=':') for x in [11.5, 15.5, 39.5, 43.5]]
    #     [ax[1].axvline(x, color='gray', ls=':') for x in [1.5, 9.5]]
    #     ax[1].set_xticklabels([0, 5, 10])
    #     ax[1].set_xlabel('time (s)')
    #     fig.savefig(f'raw_fr_idx{i}_{idx[i]}.png',bbox_inches='tight')
    #     # plt.show()


def refine_svm():
    delay = 6
    fstr = np.load("ctd.npz", allow_pickle=True)
    features_per_su = fstr["features_per_su"].tolist()
    baseline_WRS_p3_p6 = baseline_statstics(features_per_su)
    sus_trans_flag = per_second_stats.process_all()  # 33172 x 4, sust,trans,switch,unclassified
    wrs_p = baseline_WRS_p3_p6[:, 0] if delay == 3 else baseline_WRS_p3_p6[:, 1]
    trans_feat = [features_per_su[i] for i in np.nonzero(np.logical_and(sus_trans_flag[1, :], wrs_p > 0.05))[0]]
    features_per_su_all = features_per_su
    features_per_su = trans_feat
    one_cv = cross_time_decoding_actual(trans_feat, n_neuron=1000, n_trial=(20, 25), delay=6,
                                        bin_range=np.arange(4, 12))


# %% main
if __name__ == "__main__":

    (sus50, trans50, trans1000) = ctd_correct_error_all(denovo=True, to_plot=True, delay=3, cpu=30)
    (sus50, trans50, trans1000) = ctd_correct_error_all(denovo=True, to_plot=True, delay=6, cpu=30)

    cross_time_decoding(denovo=False, to_plot=True, delay=6, cpu=30)
    cross_time_decoding(denovo=False, to_plot=True, delay=3, cpu=30)

    ctd_correct_error_all(denovo=True, to_plot=True, delay=3, cpu=30)
    ctd_correct_error_all(denovo=True, to_plot=True, delay=6, cpu=30)

    sys.exit()
    

    (sus, trans) = last_bin_decoding()
    np.savez_compressed('sus_transient_decoding.npz', sus=sus, trans=trans)
    sus_per_count = []
    for count_bin in sus:
        mm = np.mean([x[0] for x in count_bin])
        std = np.std([x[0] for x in count_bin])
        mm_shuf = np.mean([x[2] for x in count_bin])
        std_shuf = np.std([x[2] for x in count_bin])
        sus_per_count.append([mm, std, mm_shuf, std_shuf])
    trans_per_count = []
    for count_bin in trans:
        mm = np.mean([x[0] for x in count_bin])
        std = np.std([x[0] for x in count_bin])
        mm_shuf = np.mean([x[2] for x in count_bin])
        std_shuf = np.std([x[2] for x in count_bin])
        trans_per_count.append([mm, std, mm_shuf, std_shuf])

    (fh, axes) = plt.subplots(1, 1, figsize=(7, 5), dpi=300)
    mmsus = list(zip(*sus_per_count))[0]
    mmts = list(zip(*trans_per_count))[0]
    mmshuf = list(zip(*trans_per_count))[2]

    stdsus = list(zip(*sus_per_count))[1]
    stdts = list(zip(*trans_per_count))[1]
    stdshuf = list(zip(*trans_per_count))[3]

    axes.fill_between(np.arange(len(mmsus)),
                      [mmsus[i] + stdsus[i] for i in np.arange(len(mmsus))],
                      [mmsus[i] - stdsus[i] for i in np.arange(len(mmsus))], fc='r', alpha=0.2)
    axes.fill_between(np.arange(len(mmts)),
                      [mmts[i] + stdts[i] for i in np.arange(len(mmts))],
                      [mmts[i] - stdts[i] for i in np.arange(len(mmts))], fc='b', alpha=0.2)
    axes.fill_between(np.arange(len(mmshuf)),
                      [mmshuf[i] + stdshuf[i] for i in np.arange(len(mmshuf))],
                      [mmshuf[i] - stdshuf[i] for i in np.arange(len(mmshuf))], fc='k', alpha=0.2)

    ph0 = axes.plot(mmsus, 'r-')
    ph1 = axes.plot(mmts, 'b-')
    ph2 = axes.plot(mmshuf, 'k-')

    axes.set_xticks(np.arange(-1, 40, 4))
    axes.set_xticklabels(np.arange(0, 1001, 100))
    axes.set_ylim((40, 100))
    axes.set_ylabel('decoding accuracy')
    axes.set_xlabel('Number of single units')
    axes.legend((ph0[0], ph1[0], ph2[0]), ('sustained (157/33172)', 'transient (8680/33172)', 'shuffled'))
    plt.show()
    fh.savefig('n_su_vs_decoding_accu.png')
