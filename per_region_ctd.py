# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:49:17 2020

@author: Libra
"""
import os

# import sys
import csv
import h5py
import numpy as np
import scipy.stats as stats
from sklearn.svm import SVC
from sklearn.model_selection import KFold
from mcc.mcc import MaximumCorrelationClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import su_region_align as align
import per_second_stats
from multiprocessing import Pool
from matplotlib import rcParams


class ctd_stats:
    def __init__(self):
        self.su_list = []
        self.use_ranksum = True

    def bool_stats_test(self, A, B):
        ### TODO alternative use of perm_test
        if self.use_ranksum:
            try:
                (_stat, p) = stats.mannwhitneyu(A.flatten(), B.flatten(), alternative="two-sided")
                return p < 0.001

            except ValueError:
                return False
        else:
            return self.exact_mc_perm_test(A, B, 1000)

    def exact_mc_perm_test(self, xs, ys, nmc):
        n, k = len(xs), 0
        diff = np.abs(np.mean(xs) - np.mean(ys))
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
        return k / nmc

    ### unprocessed data entry point

    def processCTDStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None):

        ### TODO: when variables are empty
        # [bin:trial:SU]
        # breakpoint()
        FR_D3_S1 = trial_FR[
            :, np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 4, welltrain_window, correct_resp,)), axis=0,), :,
        ]
        FR_D3_S2 = trial_FR[
            :, np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 8, welltrain_window, correct_resp,)), axis=0,), :,
        ]
        FR_D6_S1 = trial_FR[
            :, np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 4, welltrain_window, correct_resp,)), axis=0,), :,
        ]
        FR_D6_S2 = trial_FR[
            :, np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 8, welltrain_window, correct_resp,)), axis=0,), :,
        ]

        FR_D3_S1_ERR = trial_FR[
            :, np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 4, np.logical_not(correct_resp))), axis=0,), :,
        ]
        FR_D3_S2_ERR = trial_FR[
            :, np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 8, np.logical_not(correct_resp))), axis=0,), :,
        ]
        FR_D6_S1_ERR = trial_FR[
            :, np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 4, np.logical_not(correct_resp))), axis=0,), :,
        ]
        FR_D6_S2_ERR = trial_FR[
            :, np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 8, np.logical_not(correct_resp))), axis=0,), :,
        ]

        for su_idx in range(trial_FR.shape[2]):
            self.su_list.append(
                {
                    "S1_3": np.squeeze(FR_D3_S1[:, :, su_idx]),
                    "S2_3": np.squeeze(FR_D3_S2[:, :, su_idx]),
                    "S1_6": np.squeeze(FR_D6_S1[:, :, su_idx]),
                    "S2_6": np.squeeze(FR_D6_S2[:, :, su_idx]),
                    "S1_3_ERR": np.squeeze(FR_D3_S1_ERR[:, :, su_idx]),
                    "S2_3_ERR": np.squeeze(FR_D3_S2_ERR[:, :, su_idx]),
                    "S1_6_ERR": np.squeeze(FR_D6_S1_ERR[:, :, su_idx]),
                    "S2_6_ERR": np.squeeze(FR_D6_S2_ERR[:, :, su_idx]),
                }
            )

    def get_features(self):
        # breakpoint()
        return self.su_list


### Main program entry point


def get_dataset(denovo):
    features_per_su = []
    reg_list = []
    if denovo:
        dpath = align.get_root_path()
        for path in align.traverse(dpath):
            print(path)

            # SU_ids = []
            trial_FR = None
            trials = None
            if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
                continue
            done_read = False
            while not done_read:
                try:
                    with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
                        # print(list(ffr.keys()))
                        if not "SU_id" in ffr.keys():
                            done_read = True
                            print("missing su_id key in path ", path)
                            continue
                        # dset = ffr["SU_id"]
                        # SU_ids = np.array(dset, dtype="uint16")
                        dset = ffr["FR_All"]
                        trial_FR = np.array(dset, dtype="double")
                        dset = ffr["Trials"]
                        trials = np.array(dset, dtype="double").T
                    done_read = True
                except OSError:
                    print("h5py read error handled")

            if trials is None:
                continue

            (_perf_desc, perf_code, welltrain_window, correct_resp,) = align.judgePerformance(trials)

            if perf_code != 3:
                continue

            suid_reg = []
            with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
                l = list(csv.reader(csvfile))[1:]
                suid_reg = [list(i) for i in zip(*l)]

            currStats = ctd_stats()
            currStats.processCTDStats(trial_FR, trials, welltrain_window, correct_resp)
            features_per_su.extend(currStats.get_features())
            reg_list.extend(suid_reg[1])

            ### DEBUG in small subsets
            # if len(features_per_su)>50:
            #     break

        ### save to npz file
        np.savez_compressed("ctd.npz", features_per_su=features_per_su, reg_list=reg_list)

    ### load from npz file
    else:
        fstr = np.load("ctd.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        reg_list = fstr["reg_list"].tolist()

    return (features_per_su, reg_list)


def same_time_decoding(denovo, n_sel, delay=3, limit_bins=None):
    (features_per_su, _reg_list) = get_dataset(denovo)

    scaler = MinMaxScaler()
    avail_sel = [(x["S1_3"].shape[1] >= n_sel and x["S2_3"].shape[1] >= n_sel) for x in features_per_su]
    clf = SVC(kernel="linear")  # bins, trials
    one_dir = []
    rng = np.random.default_rng()
    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]

    if limit_bins is None:
        limit_bins = features_per_su[0]["S1_3"].shape[0]

    for bin_idx in range(limit_bins):
        X1 = np.vstack(tuple([su[keys[0]][bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])).T
        X2 = np.vstack(tuple([su[keys[1]][bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])).T
        y1 = np.ones((X1.shape[0]))
        y2 = np.zeros_like(y1)
        y = np.hstack((y1, y2))
        y_shuf = y.copy()
        rng.shuffle(y_shuf)
        X = scaler.fit_transform(np.vstack((X1, X2)))
        scores = cross_val_score(clf, X, y, cv=n_sel) * 100
        scores_shuffled = cross_val_score(clf, X, y_shuf, cv=n_sel) * 100
        one_dir.append(
            [np.mean(scores), np.std(scores), np.mean(scores_shuffled), np.std(scores_shuffled),]
        )

    mm = np.array(one_dir)[:, 0]
    sem = np.array(one_dir)[:, 1] / np.sqrt(scores.shape[0] - 1)
    mm_shuf = np.array(one_dir)[:, 2]
    sem_shuf = np.array(one_dir)[:, 3] / np.sqrt(scores.shape[0] - 1)
    fh = plt.figure()
    (lh,) = plt.plot(mm, color="r")
    plt.fill_between(np.arange(mm.shape[0]), mm - sem, mm + sem, color="r", alpha=0.2)
    (lhs,) = plt.plot(mm_shuf, color="k")
    plt.fill_between(
        np.arange(mm.shape[0]), mm_shuf - sem_shuf, mm_shuf + sem_shuf, color="k", alpha=0.2,
    )
    ax = plt.gca()

    if delay == 3:
        [ax.axvline(x, lw=0.5, ls=":", c="w") for x in [11.5, 15.5, 27.5, 31.5]]
    else:
        [ax.axvline(x, lw=0.5, ls=":", c="w") for x in [11.5, 15.5, 39.5, 43.5]]
    ax.set_xticks([11.5, 31.5, 51.5])
    ax.set_xticklabels([0, 5, 10])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("classification accuracy")
    ax.legend((lh, lhs), ("experiment", "shuffled"))
    fh.savefig("same_time_decoding.png", dpi=300, bbox_inches="tight")


def cross_time_decoding_old(features_per_su, n_sel, delay=3, limit_bins=None, reg_name="All"):
    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]

    scaler = MinMaxScaler()
    avail_sel = [(x[keys[0]].shape[1] >= n_sel and x[keys[1]].shape[1] >= n_sel) for x in features_per_su]

    if sum(avail_sel) < 10:
        print("Not enough SU with suffcient trials")
        return None

    clf = SVC(kernel="linear")  # bins, trials

    if limit_bins is None:
        limit_bins = features_per_su[0][keys[0]].shape[0]

    score_mat = np.zeros((limit_bins, limit_bins))
    for template_bin_idx in range(limit_bins):
        X1 = np.vstack(
            tuple([su[keys[0]][template_bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])
        ).T
        X2 = np.vstack(
            tuple([su[keys[1]][template_bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])
        ).T
        y1 = np.ones((X1.shape[0]))
        y2 = np.zeros_like(y1)
        y = np.hstack((y1, y2))
        # y_shuf=y.copy()
        # rng.shuffle(y_shuf)
        template_X = scaler.fit_transform(np.vstack((X1, X2)))
        clf.fit(template_X, y)
        for test_bin_idx in range(limit_bins):
            XX1 = np.vstack(
                tuple([su[keys[0]][test_bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])
            ).T
            XX2 = np.vstack(
                tuple([su[keys[1]][test_bin_idx, :n_sel] for (su, tf) in zip(features_per_su, avail_sel) if tf])
            ).T

            test_X = scaler.fit_transform(np.vstack((XX1, XX2)))
            score_mat[template_bin_idx, test_bin_idx] = clf.score(test_X, y)

    score_mat = score_mat * 100
    (fh, ax) = plt.subplots()
    im = plt.imshow(score_mat, cmap="jet", aspect="auto", origin="lower", vmin=50, vmax=100)
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
    fh.savefig(f"cross_time_decoding_{suffix}_{reg_name}.png", dpi=300, bbox_inches="tight")
    np.save(f"score_mat_{suffix}_{reg_name}.npy", score_mat)
    return score_mat

    ### disabled due to missing trials
    # availErrTrials=[]
    # for su in features_per_su:
    #     availErrTrials.append([su['S1_3_ERR'].shape[1],su['S2_3_ERR'].shape[1],su['S1_6_ERR'].shape[1],su['S2_6_ERR'].shape[1]])


def cross_time_decoding(
    trans_feat, to_plot=False, delay=6, repeats=10, decoder="SVC", reg="RND",
):

    # baseline_WRS_p3_p6 = baseline_statstics(features_per_su)
    # wrs_p = baseline_WRS_p3_p6[:, 0] if delay == 3 else baseline_WRS_p3_p6[:, 1]

    (avail, trans50) = ctd_actual(
        trans_feat,
        n_neuron=50,
        n_trial=(20, 25),
        delay=delay,
        bin_range=np.arange(4, 52),
        decoder=decoder,
        repeats=repeats,
    )

    if not avail:
        return (False, None)

    elif to_plot:
        (fig, ax) = plt.subplots(1, 1, figsize=[15 / 25.4, 15 / 25.4], dpi=300)

        # ax.set_ylabel("template time (s)")
        # ax[1].set_title("50 transient SU")
        im2 = ax.imshow(np.array(trans50).mean(axis=0), cmap="jet", aspect="auto", origin="lower", vmin=0, vmax=100)
        # plt.colorbar(im2, ticks=[0, 50, 100], format="%d")

        ax.set_xticks([7.5, 27.5])
        ax.set_xticklabels([])
        ax.set_xlabel(reg)
        ax.set_yticks([7.5, 27.5])
        ax.set_yticklabels([])

        if delay == 6:
            [ax.axhline(x, color="w", ls=":", lw=0.5) for x in [7.5, 11.5, 35.5, 39.5]]
            [ax.axvline(x, color="w", ls=":", lw=0.5) for x in [7.5, 11.5, 35.5, 39.5]]
        elif delay == 3:
            [ax.axhline(x, color="w", ls=":", lw=0.5) for x in [7.5, 11.5, 23.5, 27.5]]
            [ax.axvline(x, color="w", ls=":", lw=0.5) for x in [7.5, 11.5, 23.5, 27.5]]

        fig.savefig(f"ctd_{delay}_{repeats}_{reg}.pdf", bbox_inches="tight")

        plt.show()
        return (True, trans50)


def ctd_actual(features_per_su, n_neuron=50, n_trial=(20, 25), delay=6, bin_range=None, decoder="MCC", repeats=20):
    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]
    avail_sel = [(x[keys[0]].shape[1] >= n_trial[1] and x[keys[1]].shape[1] >= n_trial[1]) for x in features_per_su]

    if sum(avail_sel) < n_neuron:
        print("Not enough SU with suffcient trials")
        return (False, None)

    # bins, trials
    if bin_range is None:
        bin_range = np.arange(features_per_su[0][keys[0]].shape[0])

    scaler = MinMaxScaler()
    if decoder == "MCC":
        clf = MaximumCorrelationClassifier(n_neuron)
    elif decoder == "SVC":
        clf = SVC(kernel="linear")

    kf = KFold(10)
    one_cv = []
    for i in range(repeats):
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
                    y_test = np.hstack((np.ones_like(tests) * -1, np.ones_like(tests) * 1)).T
                    # y_shuf=y.copy()
                    # rng.shuffle(y_shuf)
                    X_test = scaler.transform(X_test)
                    if decoder == "MCC":
                        score_mat[template_bin_idx, test_bin_idx] = np.mean(
                            np.hstack(
                                (
                                    clf.predict(X1[tests, :, test_bin_idx]) <= 0,
                                    clf.predict(X2[tests, :, test_bin_idx]) > 0,
                                )
                            )
                        )
                    else:
                        score_mat[template_bin_idx, test_bin_idx] = clf.score(X_test, y_test)

            score_mat = score_mat * 100
            one_cv.append(score_mat)

    return (True, one_cv)


def ctd_par(denovo=False, delay=6, proc_n=2, repeats=2):
    rcParams["pdf.fonttype"] = 42
    rcParams["ps.fonttype"] = 42
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = ["Arial"]  # breakpoint()
    # be noted n_sel was not used until later stage, e.g. saved dataset is complete
    (features_per_su, reg_list) = get_dataset(denovo)
    trans_fstr = np.load(f"sus_trans_pie_{delay}.npz")
    # list(trans6.keys())
    # sust = trans_fstr["sust"]
    trans = trans_fstr["transient"]
    # reg_arr=trans_fstr['reg_arr']  tested to be identical to ctd.npz source

    trans50_all = []
    reg_set = sorted(set(reg_list))
    # reg_count = [reg_list.count(x) for x in reg_set]

    reg_count = [[reg, np.sum(np.logical_and([x == reg for x in reg_list], trans))] for reg in reg_set]

    reg_set = [x[0] for x in reg_count if x[1] >= 50]

    curr_pool = Pool(processes=proc_n)

    all_proc = []
    # reg_offset = (
    #     (reg_onset + proc_n) if (reg_onset <= len(reg_set) - proc_n) else len(reg_set)
    # )

    for reg_idx in range(len(reg_set)):
        curr_feat_su = [f for (f, reg, t) in zip(features_per_su, reg_list, trans) if (reg == reg_set[reg_idx] and t)]
        # trans_feat, to_plot=False, delay=6, repeats=10, decoder="SVC", reg="RND",
        all_proc.append(
            curr_pool.apply_async(
                cross_time_decoding,
                args=(curr_feat_su,),
                kwds={"to_plot": True, "delay": delay, "repeats": repeats, "decoder": "SVC", "reg": reg_set[reg_idx],},
            )
        )

        # print("set")

    for one_proc in all_proc:
        (avail, trans50) = one_proc.get()
        trans50_all.append(trans50)
        # print("get")

    curr_pool.close()
    curr_pool.join()

    trans50Avail = [x for x in trans50 if x is not None]
    reg_setAvail = [reg_set[i] for i in range(len(trans50)) if trans50[i] is not None]
    np.savez_compressed(f"per_region_ctd_{delay}_{repeats}.npz", trans50=trans50Avail, reg_set=reg_setAvail)

    return [reg_set, trans50_all]


if __name__ == "__main__":

    (reg_set, trans50) = ctd_par(denovo=True, delay=6, proc_n=8, repeats=2)
    (reg_set, trans50) = ctd_par(denovo=False, delay=6, proc_n=8, repeats=100)
