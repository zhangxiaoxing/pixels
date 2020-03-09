import per_second_stats
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt


def decoding():
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


def same_time_decoding(features_per_su, n_neuron=None, n_trial=None, delay=6, bin_range=None, repeat=10):
    if bin_range is None:
        bin_range = np.arange(36, 40)

    if n_neuron is None:
        n_neuron = 300

    if n_trial is None:
        n_trial = (20, 25)

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


def cross_time_decoding(features_per_su, n_neuron=None, n_trial=None, delay=6, bin_range=None, repeat=10):
    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]

    scaler = MinMaxScaler()
    avail_sel = [
        (x[keys[0]].shape[1] >= n_trial[1] and x[keys[1]].shape[1] >= n_trial[1])
        for x in features_per_su
    ]

    if sum(avail_sel) < 10:
        print('Not enough SU with suffcient trials')
        return None

    clf = LinearSVC()
    # bins, trials

    if limit_bins is None:
        limit_bins = features_per_su[0][keys[0]].shape[0]

    score_mat = np.zeros((limit_bins, limit_bins))
    for template_bin_idx in range(limit_bins):
        X1 = np.vstack(
            tuple([su[keys[0]][template_bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])).T
        X2 = np.vstack(
            tuple([su[keys[1]][template_bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])).T
        y1 = np.ones((X1.shape[0]))
        y2 = np.zeros_like(y1)
        y = np.hstack((y1, y2))
        # y_shuf=y.copy()
        # rng.shuffle(y_shuf)
        template_X = scaler.fit_transform(np.vstack((X1, X2)))
        clf.fit(template_X, y)
        for test_bin_idx in range(limit_bins):
            XX1 = np.vstack(
                tuple([su[keys[0]][test_bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])).T
            XX2 = np.vstack(
                tuple([su[keys[1]][test_bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])).T

            test_X = scaler.fit_transform(np.vstack((XX1, XX2)))
            score_mat[template_bin_idx, test_bin_idx] = clf.score(test_X, y)

    score_mat = score_mat * 100
    (fh, ax) = plt.subplots()
    im = plt.imshow(
        score_mat, cmap="jet", aspect="auto", origin="lower", vmin=50, vmax=100
    )
    plt.colorbar(im, ticks=[50, 75, 100], format="%d")
    suffix = None
    if delay == 3:
        [ax.axvline(x, lw=0.5, ls=":", c="w") for x in [11.5, 15.5, 27.5, 31.5]]
        [ax.axhline(x, lw=0.5, ls=":", c="w") for x in [11.5, 15.5, 27.5, 31.5]]
        suffix = "3S"
    else:
        [ax.axvline(x, lw=0.5, ls=":", c="w") for x in [11.5, 15.5, 39.5, 43.5]]
        [ax.axhline(x, lw=0.5, ls=":", c="w") for x in [11.5, 15.5, 39.5, 43.5]]
        suffix = "6S"

    ax.set_xticks([11.5, 31.5, 51.5])
    ax.set_xticklabels([0, 5, 10])
    ax.set_xlabel("Time (s)")
    ax.set_yticks([11.5, 31.5, 51.5])
    ax.set_yticklabels([0, 5, 10])
    ax.set_ylabel("Time (s)")
    ax.set_title(f"{reg_name} CTD {suffix} delay")
    fh.savefig(
        f"cross_time_decoding_{suffix}_{reg_name}.png", dpi=300, bbox_inches="tight"
    )
    np.save(f"score_mat_{suffix}_{reg_name}.npy", score_mat)
    return score_mat

    ### disabled due to missing trials
    # availErrTrials=[]
    # for su in features_per_su:
    #     availErrTrials.append([su['S1_3_ERR'].shape[1],su['S2_3_ERR'].shape[1],su['S1_6_ERR'].shape[1],su['S2_6_ERR'].shape[1]])


if __name__ == "__main__":
    (sus, trans) = decoding()
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
