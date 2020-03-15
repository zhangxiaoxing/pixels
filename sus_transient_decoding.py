import sys
import per_second_stats
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler
from multiprocessing import Pool
import matplotlib.pyplot as plt


def last_bin_decoding():
    fstr = np.load("ctd.npz", allow_pickle=True)
    features_per_su = fstr["features_per_su"].tolist()
    # reg_list = fstr["reg_list"].tolist()
    sus_trans_flag = per_second_stats.process_all()  # 33172 x 4, sust,trans,switch,unclassified
    sus_feat = [features_per_su[i] for i in np.nonzero(sus_trans_flag[:, 0])[0]]
    trans_feat = [features_per_su[i] for i in np.nonzero(sus_trans_flag[:, 1])[0]]
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
        for one_su in su_selected_features:
            trial_to_select1 = np.random.choice(one_su[keys[0]].shape[1], n_trial[0], replace=False)
            X1.append([np.mean(one_su[keys[0]][bin_range, t]) for t in trial_to_select1])
            trial_to_select2 = np.random.choice(one_su[keys[1]].shape[1], n_trial[0], replace=False)
            X2.append([np.mean(one_su[keys[1]][bin_range, t]) for t in trial_to_select2])

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


def cross_time_decoding(denovo=False, to_plot=False, delay=6):
    repeats = 1000
    if denovo:
        fstr = np.load("ctd.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        # reg_list = fstr["reg_list"].tolist()
        sus_trans_flag = per_second_stats.process_all()  # 33172 x 4, sust,trans,switch,unclassified
        sus_feat = [features_per_su[i] for i in np.nonzero(sus_trans_flag[:, 0])[0]]
        trans_feat = [features_per_su[i] for i in np.nonzero(sus_trans_flag[:, 1])[0]]
        curr_pool = Pool(processes=24)
        sus100 = []
        trans100 = []
        trans1000 = []
        sust_proc = []
        trans100_proc = []
        trans1000_proc = []
        for i in range(repeats):
            sust_proc.append(curr_pool.apply_async(cross_time_decoding_calc, args=(sus_feat,),
                                                   kwds={"n_neuron": 100, "n_trial": (20, 25), "delay": delay,
                                                         "bin_range": np.arange(4, 52)}))

            trans100_proc.append(curr_pool.apply_async(cross_time_decoding_calc, args=(trans_feat,),
                                                       kwds={"n_neuron": 100, "n_trial": (20, 25), "delay": delay,
                                                             "bin_range": np.arange(4, 52)}))
            trans1000_proc.append(curr_pool.apply_async(cross_time_decoding_calc, args=(trans_feat,),
                                                        kwds={"n_neuron": 1000, "n_trial": (20, 25), "delay": delay,
                                                              "bin_range": np.arange(4, 52)}))
        for one_proc in sust_proc:
            sus100.append(one_proc.get())

        for one_proc in trans100_proc:
            trans100.append(one_proc.get())

        for one_proc in trans1000_proc:
            trans1000.append(one_proc.get())

        curr_pool.close()
        curr_pool.join()
        np.savez_compressed(f'sus_trans_ctd_{delay}_{repeats}.npz', sus100=sus100, trans100=trans100,
                            trans1000=trans1000)
        return (sus100, trans100, trans1000)


def ctd_cross_all(denovo=False, to_plot=False, delay=3, cpu=30):
    repeats = 100
    if denovo:
        fstr = np.load("ctd.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        # reg_list = fstr["reg_list"].tolist()
        sus_trans_flag = per_second_stats.process_all(denovo=True,delay=delay).T  # 33172 x 4, sust,trans,switch,unclassified
        sus_feat = [features_per_su[i] for i in np.nonzero(sus_trans_flag[:, 0])[0]]
        trans_feat = [features_per_su[i] for i in np.nonzero(sus_trans_flag[:, 1])[0]]
        curr_pool = Pool(processes=cpu)
        sus100 = []
        trans100 = []
        trans1000 = []
        sust_proc = []
        trans100_proc = []
        trans1000_proc = []
        for i in range(repeats):
            sust_proc.append(curr_pool.apply_async(cross_time_decoding_cross, args=(sus_feat,),
                                                   kwds={"n_neuron": 100, "n_trial": (20, 25), "template_delay": delay,
                                                         "bin_range": np.arange(4, 52)}))

            trans100_proc.append(curr_pool.apply_async(cross_time_decoding_cross, args=(trans_feat,),
                                                       kwds={"n_neuron": 100, "n_trial": (20, 25), "template_delay": delay,
                                                             "bin_range": np.arange(4, 52)}))
            trans1000_proc.append(curr_pool.apply_async(cross_time_decoding_cross, args=(trans_feat,),
                                                        kwds={"n_neuron": 1000, "n_trial": (20, 25), "template_delay": delay,
                                                              "bin_range": np.arange(4, 52)}))
        for one_proc in sust_proc:
            sus100.append(one_proc.get())

        for one_proc in trans100_proc:
            trans100.append(one_proc.get())

        for one_proc in trans1000_proc:
            trans1000.append(one_proc.get())

        curr_pool.close()
        curr_pool.join()
        np.savez_compressed(f'sus_trans_ctd_{delay}_{repeats}.npz', sus100=sus100, trans100=trans100,
                            trans1000=trans1000)
    else:
        fstr = np.load(f"sus_trans_ctd_{delay}_{repeats}.npz")
        sus100 = fstr['sus100']
        trans100 = fstr['trans100']
        trans1000 = fstr['trans1000']

    sus_mat = np.array(sus100).mean(axis=0).mean(axis=0)
    trans100_mat = np.array(trans100).mean(axis=0).mean(axis=0)
    trans1000_mat = np.array(trans1000).mean(axis=0).mean(axis=0)
    if to_plot:
        (fig, ax) = plt.subplots(1, 3, figsize=[9, 3], dpi=200)
        ax[0].imshow(sus_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[0].set_title('100 sustained SU')
        ax[0].set_ylabel('scoring time (s)')
        ax[1].imshow(trans100_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[1].set_title('100 transient SU')
        im2 = ax[2].imshow(trans1000_mat, cmap="jet", aspect="auto", origin='lower', vmin=0, vmax=100)
        ax[2].set_title('1000 transient SU')
        plt.colorbar(im2, ticks=[0, 50, 100], format="%d")

        for oneax in ax:
            oneax.set_xticks([7.5, 27.5])
            oneax.set_xticklabels([0, 5])
            oneax.set_xlabel('template time (s)')
            oneax.set_yticks([7.5, 27.5])
            oneax.set_yticklabels([0, 5])

            if delay == 6:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
            elif delay == 3:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]

        fig.savefig('ctd_cross_delay_dur.png', bbox_inches='tight')
        plt.show()

    return (sus100, trans100, trans1000)


def cross_time_decoding_calc(features_per_su, n_neuron=300, n_trial=(20, 25), delay=6, bin_range=None):
    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]
    avail_sel = [(x[keys[0]].shape[1] >= n_trial[1] and x[keys[1]].shape[1] >= n_trial[1]) for x in features_per_su]

    if sum(avail_sel) < n_neuron:
        print('Not enough SU with suffcient trials')
        return None

    # bins, trials
    if bin_range is None:
        bin_range = np.arange(features_per_su[0][keys[0]].shape[0])

    scaler = MinMaxScaler()
    clf = LinearSVC()
    kf = KFold(10)
    su_index = np.random.choice(np.nonzero(avail_sel)[0], n_neuron, replace=False)
    su_selected_features = [features_per_su[i] for i in su_index]
    X1 = []
    X2 = []
    for one_su in su_selected_features:
        trial_to_select1 = np.random.choice(one_su[keys[0]].shape[1], n_trial[0], replace=False)
        X1.append([one_su[keys[0]][bin_range, t] for t in trial_to_select1])
        trial_to_select2 = np.random.choice(one_su[keys[1]].shape[1], n_trial[0], replace=False)
        X2.append([one_su[keys[1]][bin_range, t] for t in trial_to_select2])

    X1 = np.array(X1).transpose((1, 0, 2))
    X2 = np.array(X2).transpose((1, 0, 2))  # trial, SU, bin

    y1 = np.ones((X1.shape[0]))
    y2 = np.zeros_like(y1)
    y = np.hstack((y1, y2))
    one_cv = []
    for (templates, tests) in kf.split(X1):
        score_mat = np.zeros((bin_range.shape[0], bin_range.shape[0]))
        for template_bin_idx in np.arange(bin_range.shape[0]):
            X_templates = np.vstack((X1[templates, :, template_bin_idx], X2[templates, :, template_bin_idx]))
            y_templates = np.hstack((np.zeros_like(templates), np.ones_like(templates))).T
            scaler = scaler.fit(X_templates)
            X_templates = scaler.transform(X_templates)
            for test_bin_idx in np.arange(bin_range.shape[0]):
                X_test = np.vstack((X1[tests, :, test_bin_idx], X2[tests, :, test_bin_idx]))
                y_test = np.hstack((np.zeros_like(tests), np.ones_like(tests))).T
                # y_shuf=y.copy()
                # rng.shuffle(y_shuf)
                clf.fit(X_templates, y_templates)
                X_test = scaler.transform(X_test)
                score_mat[template_bin_idx, test_bin_idx] = clf.score(X_test, y_test)

        score_mat = score_mat * 100
        one_cv.append(score_mat)

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
    for one_su in su_selected_features:
        template_trial_to_select1 = np.random.choice(one_su[template_keys[0]].shape[1], n_trial[0], replace=False)
        score_trial_to_select1 = np.random.choice(one_su[score_keys[0]].shape[1], n_trial[0], replace=False)
        template_X1.append([one_su[template_keys[0]][bin_range, t] for t in template_trial_to_select1])
        score_X1.append([one_su[score_keys[0]][bin_range, t] for t in score_trial_to_select1])
        template_trial_to_select2 = np.random.choice(one_su[template_keys[1]].shape[1], n_trial[0], replace=False)
        score_trial_to_select2 = np.random.choice(one_su[score_keys[1]].shape[1], n_trial[0], replace=False)
        template_X2.append([one_su[template_keys[1]][bin_range, t] for t in template_trial_to_select2])
        score_X2.append([one_su[score_keys[1]][bin_range, t] for t in score_trial_to_select2])

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


if __name__ == "__main__":
    (sus100, trans100, trans1000) = ctd_cross_all(denovo=True,to_plot=True, delay=3, cpu=30)

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
