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
import csv
from sus_transient_decoding import ctd_actual, ctd_correct_error


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
                   :,
                   np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 4, welltrain_window, correct_resp,)), axis=0, ),
                   :,
                   ]
        FR_D3_S2 = trial_FR[
                   :,
                   np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 8, welltrain_window, correct_resp,)), axis=0, ),
                   :,
                   ]
        FR_D6_S1 = trial_FR[
                   :,
                   np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 4, welltrain_window, correct_resp,)), axis=0, ),
                   :,
                   ]
        FR_D6_S2 = trial_FR[
                   :,
                   np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 8, welltrain_window, correct_resp,)), axis=0, ),
                   :,
                   ]

        FR_D3_S1_ERR = trial_FR[
                       :, np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 4, np.logical_not(correct_resp))),
                                 axis=0, ), :,
                       ]
        FR_D3_S2_ERR = trial_FR[
                       :, np.all(np.vstack((trials[:, 5] == 3, trials[:, 2] == 8, np.logical_not(correct_resp))),
                                 axis=0, ), :,
                       ]
        FR_D6_S1_ERR = trial_FR[
                       :, np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 4, np.logical_not(correct_resp))),
                                 axis=0, ), :,
                       ]
        FR_D6_S2_ERR = trial_FR[
                       :, np.all(np.vstack((trials[:, 5] == 6, trials[:, 2] == 8, np.logical_not(correct_resp))),
                                 axis=0, ), :,
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
            [np.mean(scores), np.std(scores), np.mean(scores_shuffled), np.std(scores_shuffled), ]
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


def cross_time_decoding(
        trans_feat, to_plot=False, n_neuron=50, delay=6, repeats=10, decoder="SVC", reg="RND", correct_error=None
):
    # baseline_WRS_p3_p6 = baseline_statstics(features_per_su)
    # wrs_p = baseline_WRS_p3_p6[:, 0] if delay == 3 else baseline_WRS_p3_p6[:, 1]
    if correct_error is None:
        (avail, decode_mat, actual_number) = ctd_actual(
            trans_feat,
            n_neuron=n_neuron,
            n_trial=(20, 25),
            delay=delay,
            bin_range=np.arange(8, 40),
            decoder=decoder,
            repeats=repeats,
        )
    else:
        (avail, decode_mat, actual_number) = ctd_correct_error(
            trans_feat,
            n_neuron=n_neuron,
            n_trial=(20, 25, 2, 4),
            template_delay=6,
            bin_range=np.arange(8, 40),
            decoder=decoder,
            repeats=repeats,
        )

    if not avail:
        return (False, reg, None, None)

    elif to_plot:
        (fig, ax) = plt.subplots(1, 1, figsize=[16 / 25.4, 16 / 25.4], dpi=300)

        # ax.set_ylabel("template time (s)")
        # ax[1].set_title("50 transient SU")
        im2 = ax.imshow(np.array(decode_mat).mean(axis=0), cmap="jet", aspect="auto", origin="lower", vmin=25, vmax=75)
        # plt.colorbar(im2, ticks=[0, 50, 100], format="%d")
        ax.tick_params(direction='in', length=2, width=0.5, colors='k')
        ax.set_xticks([3.5, 23.5])
        ax.set_xticklabels([])
        ax.set_xlabel(f"{reg}_{actual_number}SU")
        ax.set_yticks([3.5, 23.5])
        ax.set_yticklabels([])

        if delay == 6:
            [ax.axhline(x, color="k", ls=":", lw=0.5) for x in [3.5, 7.5]]
            [ax.axvline(x, color="k", ls=":", lw=0.5) for x in [3.5, 7.5]]
        elif delay == 3:
            [ax.axhline(x, color="k", ls=":", lw=0.5) for x in [3.5, 7.5]]
            [ax.axvline(x, color="k", ls=":", lw=0.5) for x in [3.5, 7.5]]

        if correct_error is None:
            suffix=''
        else:
            suffix='_error'

        fig.savefig(f"ctd_{suffix}{delay}_{repeats}_{reg}_{decoder}.pdf", bbox_inches="tight")
        fig.savefig(f"ctd_{suffix}{delay}_{repeats}_{reg}_{decoder}.png", bbox_inches="tight")

        if rcParams['backend'] != 'pdf':
            plt.show()
    return (True, reg, decode_mat, actual_number)


def ctd_par(denovo=False, delay=6, proc_n=2, repeats=2, only_trans=False, n_neuron=50, decoder="SVC",
            correct_error=None):
    # be noted n_sel was not used until later stage, e.g. saved dataset is complete
    (features_per_su, reg_list) = get_dataset(denovo)
    trans_fstr = np.load(f"sus_trans_pie_{delay}.npz")
    # list(trans6.keys())
    # sust = trans_fstr["sust"]
    trans = trans_fstr["transient"]
    # reg_arr=trans_fstr['reg_arr']  tested to be identical to ctd.npz source

    decode_mat_all = []
    actual_number_all = []
    reg_set = sorted(set(reg_list))
    # reg_count = [reg_list.count(x) for x in reg_set]

    reg_count = [[reg, np.sum(np.logical_and([x == reg for x in reg_list], trans))] for reg in reg_set]

    reg_set = [x[0] for x in reg_count if x[1] >= 50]
    if proc_n > 1:
        curr_pool = Pool(processes=proc_n)

    all_proc = []
    # reg_offset = (
    #     (reg_onset + proc_n) if (reg_onset <= len(reg_set) - proc_n) else len(reg_set)
    # )

    for reg_idx in range(len(reg_set)):

        if only_trans:
            curr_feat_su = [f for (f, reg, t) in zip(features_per_su, reg_list, trans) if
                            (reg == reg_set[reg_idx] and t)]
        else:
            curr_feat_su = [f for (f, reg) in zip(features_per_su, reg_list) if reg == reg_set[reg_idx]]

        # trans_feat, to_plot=False, delay=6, repeats=10, decoder="SVC", reg="RND",
        if proc_n > 1:
            all_proc.append(
                curr_pool.apply_async(
                    cross_time_decoding,
                    args=(curr_feat_su,),
                    kwds={"to_plot": True, "n_neuron": n_neuron, "delay": delay, "repeats": repeats, "decoder": decoder,
                          "reg": reg_set[reg_idx], "correct_error": ctd_correct_error, },
                )
            )
        else:
            (avail, reg, decode_mat, actual_number) = cross_time_decoding(curr_feat_su,
                                                                          to_plot=True, n_neuron=n_neuron, delay=delay,
                                                                          repeats=repeats,
                                                                          reg=reg_set[reg_idx],
                                                                          correct_error=ctd_correct_error,
                                                                          )
            decode_mat_all.append(decode_mat)
            actual_number_all.append(actual_number)

        # print("set")
    if proc_n > 1:
        for one_proc in all_proc:
            (avail, reg, decode_mat, actual_number) = one_proc.get()
            decode_mat_all.append(decode_mat)
            actual_number_all.append(actual_number)
            # print("get")
        curr_pool.close()
        curr_pool.join()

    decode_avail = [x for x in decode_mat_all if x is not None]
    reg_set_avail = [reg_set[i] for i in range(len(decode_mat_all)) if decode_mat_all[i] is not None]
    actual_number_avail = [actual_number_all[i] for i in range(len(decode_mat_all)) if decode_mat_all[i] is not None]
    if correct_error is None:
        suffix = ''
    else:
        suffix = '_error'
    np.savez_compressed(f"per_region_ctd_{delay}_{repeats}_{decoder}{suffix}.npz", decode_mat=decode_avail,
                        reg_set=reg_set_avail,
                        actual_number=actual_number_avail)

    return [reg_set, decode_mat_all]


def per_region_corr():
    fstr = np.load('per_region_ctd_6_100_SVC.npz')
    # list(fstr.keys())
    # Out[3]: ['trans50', 'reg_set']
    decode_mat = fstr['decode_mat']
    actual_number = fstr['actual_number']
    reg_set = fstr['reg_set']
    reg_list = reg_set.tolist()

    exclude = [reg_list.index('Unlabeled'), reg_list.index('int')]
    decode_mat = decode_mat[np.hstack((np.arange(23), 24)), :, :, :]

    reg_dist = []
    # early, late, diag, off-diag
    for i in range(decode_mat.shape[0]):
        early = np.mean([decode_mat[i, :, x, x] for x in range(8, 20)])
        late = np.mean([decode_mat[i, :, x, x] for x in range(20, 32)])
        diag = np.mean([early, late])
        offdiagMat = np.mean(np.squeeze(decode_mat[i, :, 8:32, 8:32]), axis=0)
        for odi in range(offdiagMat.shape[0]):
            offdiagMat[odi, odi] = np.nan
        offdiag = np.nanmean(offdiagMat)
        reg_dist.append([early, late, diag, offdiag])

    reg_dist_arr = np.array(reg_dist)
    for tidx in range(decode_mat.shape[0]):
        reg_dist[tidx].insert(0, actual_number[tidx])
        reg_dist[tidx].insert(0, reg_list[tidx])

    with open("all_neuron_decoding.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in reg_dist:
            cwriter.writerow(row)

    (fig, ax) = plt.subplots(1, 1, figsize=(17 / 2.54, 8 / 2.54), dpi=300)

    for tidx in range(decode_mat.shape[0]):
        ax.plot(reg_dist_arr[tidx, 0], reg_dist_arr[tidx, 1], 'r.')
        ax.text(reg_dist_arr[tidx, 0], reg_dist_arr[tidx, 1], f'{reg_list[tidx]}')
        # plt.text(tidx,tidx,'test')
    ax.plot([0, 100], [0, 100], '--k', lw=0.5)
    ax.set_xlim([50, 90])
    ax.set_ylim([50, 80])
    ax.set_xlabel('early delay decoding accuracy')
    ax.set_ylabel('late delay')
    ax.set_xticks([50, 60, 70, 80, 90])
    ax.set_yticks([50, 60, 70, 80])
    fig.savefig('per_region_early_late.pdf', bbox_inches='tight')
    if rcParams['backend'] != 'pdf':
        plt.show()

    (fig, ax) = plt.subplots(1, 1, figsize=(17 / 2.54, 8 / 2.54), dpi=300)

    for tidx in range(decode_mat.shape[0]):
        ax.plot(reg_dist_arr[tidx, 2], reg_dist_arr[tidx, 3], 'r.')
        ax.text(reg_dist_arr[tidx, 2], reg_dist_arr[tidx, 3], f'{reg_list[tidx]}')
        # plt.text(tidx,tidx,'test')
    ax.plot([0, 100], [0, 100], '--k', lw=0.5)
    ax.set_xlim([50, 80])
    ax.set_ylim([50, 80])
    ax.set_xlabel('diagonal decoding accuracy')
    ax.set_ylabel('off-diagonal')
    ax.set_xticks([50, 60, 70, 80])
    ax.set_yticks([50, 60, 70, 80])
    fig.savefig('per_region_diag_off.pdf', bbox_inches='tight')
    if rcParams['backend'] != 'pdf':
        plt.show()


if __name__ == "__main__":
    # rcParams['backend'] = 'pdf'
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    # (reg_set, trans50) = ctd_par(denovo=False, delay=6, proc_n=8, repeats=2)
    (reg_set, all_su) = ctd_par(denovo=False, delay=6, proc_n=64, repeats=100, only_trans=False, n_neuron=None,
                                decoder="SVC", correct_error='error')
    # per_region_corr()
    # per_region_corr()
