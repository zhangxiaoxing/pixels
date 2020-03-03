# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 12:03:06 2020

@author: Libra
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as stats
import sklearn.metrics as metrics


class zxStats:
    def __init__(self):
        self.regionName = ""
        self.heatmapL = []
        self.heatmapR = []

        self.test_heatmapL = []
        self.test_heatmapR = []

        self.selectivity = []

        self.statistical_sample_selective = []

        self.statistical_test_selective = []

        self.statistical_pair_selective = []

        self.statistical_hitMiss_selective = []
        self.statistical_CRFalse_selective = []

        self.modulation = []

        self.auroc = []

        self.pair_heatmapL = []
        self.pair_heatmapR = []
        self.pair_selectivity = []

        self.featureVector = []

        self.row_sel_6 = np.concatenate(
            (np.arange(16), np.arange(16, 40, 2), np.arange(40, 68))
        )

        self.row_sel_3 = np.arange(56)

        self.use_ranksum = True

    def gauss_average(self, x):
        return np.convolve(x, [0.1968, 0.6063, 0.1968], "same")

    def baselineVector(self, one_su_trial_FR):
        base = one_su_trial_FR[:, 0:8].flatten()
        if np.std(base):
            return (np.mean(base), np.std(base))
        print("Constant baseline")
        return (np.mean(base), 0.5)

    def addTrialFRs(self, trial_FR, trials, su_sel=[], welltrain_window=[], correctResp=[]):
        if isinstance(su_sel, np.ndarray):
            trial_FR_sel = trial_FR[:, :, su_sel]
        else:
            trial_FR_sel = trial_FR

        perf_sel = np.logical_and(correctResp, welltrain_window)
        trial_sel = trials[perf_sel, :]
        trial_FR_perf_sel = trial_FR_sel[:, perf_sel, :]

        self.addSampleSelect(trial_FR_perf_sel, trial_sel)

        self.addPairSelect(trial_FR_perf_sel, trial_sel)

        self.addStatisticalSampleSel(trial_FR_perf_sel, trial_sel)

        self.addStatisticalTestSel(trial_FR_perf_sel, trial_sel)

        self.addHitMissSel(trial_FR, trials, su_sel, perf_sel)
        self.addCRFalseSel(trial_FR, trials, su_sel, perf_sel)

    def get_normalized_fr(self, trial_FR, trials, welltrain_window=None, correctResp=None):

        ### TODO: get normalized FR for epys scaling
        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(trial_FR[:, (welltrain_window
                                            & correctResp
                                            & (trials[:, 5] == 6)), su_idx]).T
            (base_mean, base_std) = self.baselineVector(onesu)
            allbins = np.mean((onesu - base_mean) / base_std, axis=0)
            self.modulation.append(
                np.abs([np.mean(allbins[i:i + 4]) for i in np.arange(4, 56, 4)])
            )
        return np.array(self.modulation).T

    def addHitMissSel(self, trial_FR, trials, su_sel, perf_sel):
        trial_FR_sel = trial_FR[:, :, su_sel]

        trial_sel_hit_3 = np.logical_and(
            perf_sel, np.logical_and(trials[:, 2] != trials[:, 3], trials[:, 5] == 3)
        )
        trial_sel_hit_6 = np.logical_and(
            perf_sel, np.logical_and(trials[:, 2] != trials[:, 3], trials[:, 5] == 6)
        )

        # Miss
        trial_sel_miss_3 = np.logical_and(
            np.logical_and(trials[:, 2] != trials[:, 3], trials[:, 4] < 0),
            trials[:, 5] == 3,
        )

        trial_sel_miss_6 = np.logical_and(
            np.logical_and(trials[:, 2] != trials[:, 3], trials[:, 4] < 0),
            trials[:, 5] == 6,
        )

        for su_idx in range(trial_FR_sel.shape[2]):
            onesu = np.squeeze(trial_FR_sel[:, :, su_idx]).T
            hit_trials_3 = onesu[trial_sel_hit_3, :][:, self.row_sel_3]
            miss_trials_3 = onesu[trial_sel_miss_3, :][:, self.row_sel_3]
            # scaled to 3s delay, e.g. 56 bins
            hit_trials_6 = onesu[trial_sel_hit_6, :][:, self.row_sel_6]

            miss_trials_6 = onesu[trial_sel_miss_6, :][:, self.row_sel_6]

            left_trials = np.concatenate((hit_trials_3, hit_trials_6))
            right_trials = np.concatenate((miss_trials_3, miss_trials_6))

            bins = np.ones(left_trials.shape[1])
            for bin_idx in range(left_trials.shape[1]):
                try:
                    (stat, p) = stats.mannwhitneyu(
                        left_trials[:, bin_idx].flatten(),
                        right_trials[:, bin_idx].flatten(),
                        alternative="two-sided",
                    )
                    bins[bin_idx] = p < 0.05
                except ValueError:
                    bins[bin_idx] = 0
            self.statistical_hitMiss_selective.append(bins)

    def addCRFalseSel(self, trial_FR, trials, su_sel, perf_sel):
        trial_FR_sel = trial_FR[:, :, su_sel]

        trial_sel_CR_3 = np.logical_and(
            perf_sel, np.logical_and(trials[:, 2] == trials[:, 3], trials[:, 5] == 3)
        )
        trial_sel_CR_6 = np.logical_and(
            perf_sel, np.logical_and(trials[:, 2] == trials[:, 3], trials[:, 5] == 6)
        )

        # Miss
        trial_sel_false_3 = np.logical_and(
            np.logical_and(trials[:, 2] == trials[:, 3], trials[:, 4] > 0),
            trials[:, 5] == 3,
        )

        trial_sel_false_6 = np.logical_and(
            np.logical_and(trials[:, 2] == trials[:, 3], trials[:, 4] > 0),
            trials[:, 5] == 6,
        )
        if trial_sel_false_3.sum() + trial_sel_false_6.sum() < 10:
            return

        for su_idx in range(trial_FR_sel.shape[2]):
            onesu = np.squeeze(trial_FR_sel[:, :, su_idx]).T
            CR_trials_3 = onesu[trial_sel_CR_3, :][:, self.row_sel_3]
            false_trials_3 = onesu[trial_sel_false_3, :][:, self.row_sel_3]
            # scaled to 3s delay, e.g. 56 bins
            CR_trials_6 = onesu[trial_sel_CR_6, :][:, self.row_sel_6]

            false_trials_6 = onesu[trial_sel_false_6, :][:, self.row_sel_6]

            left_trials = np.concatenate((CR_trials_3, CR_trials_6))
            right_trials = np.concatenate((false_trials_3, false_trials_6))

            bins = np.ones(left_trials.shape[1])
            for bin_idx in range(left_trials.shape[1]):
                try:
                    (stat, p) = stats.mannwhitneyu(
                        left_trials[:, bin_idx].flatten(),
                        right_trials[:, bin_idx].flatten(),
                        alternative="two-sided",
                    )
                    bins[bin_idx] = p < 0.05
                except ValueError:
                    bins[bin_idx] = 0
            self.statistical_CRFalse_selective.append(bins)

    def addStatisticalSampleSel(self, trial_FR, trial_sel):
        trial_sel_left_3 = np.logical_and(trial_sel[:, 2] == 4, trial_sel[:, 5] == 3)
        trial_sel_left_6 = np.logical_and(trial_sel[:, 2] == 4, trial_sel[:, 5] == 6)
        trial_sel_right_3 = np.logical_and(trial_sel[:, 2] == 8, trial_sel[:, 5] == 3)
        trial_sel_right_6 = np.logical_and(trial_sel[:, 2] == 8, trial_sel[:, 5] == 6)

        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(trial_FR[:, :, su_idx]).T

            left_trials_3 = onesu[trial_sel_left_3, :][:, self.row_sel_3]
            right_trials_3 = onesu[trial_sel_right_3, :][:, self.row_sel_3]

            left_trials_6 = onesu[trial_sel_left_6, :][:, self.row_sel_6]
            right_trials_6 = onesu[trial_sel_right_6, :][:, self.row_sel_6]

            left_trials = np.concatenate((left_trials_3, left_trials_6))

            right_trials = np.concatenate((right_trials_3, right_trials_6))

            bins = np.ones(56)

            auc = np.zeros(56)

            for bin_idx in range(left_trials.shape[1]):
                try:
                    (stat, p) = stats.mannwhitneyu(
                        left_trials[:, bin_idx].flatten(),
                        right_trials[:, bin_idx].flatten(),
                        alternative="two-sided",
                    )
                    bins[bin_idx] = p < 0.05

                except ValueError:
                    bins[bin_idx] = 0

                try:
                    auc[bin_idx] = metrics.roc_auc_score(
                        np.concatenate((
                            np.zeros(left_trials.shape[0]),
                            np.ones(right_trials.shape[0])
                        )),
                        np.concatenate((
                            left_trials[:, bin_idx].flatten(),
                            right_trials[:, bin_idx].flatten(),
                        ))
                    )
                except ValueError:
                    auc[bin_idx] = 0.5

            self.statistical_sample_selective.append(bins)

            self.auroc.append(self.gauss_average(auc))

    def addStatisticalTestSel(self, trial_FR, trial_sel):
        trial_sel_left_3 = np.logical_and(trial_sel[:, 3] == 4, trial_sel[:, 5] == 3)
        trial_sel_left_6 = np.logical_and(trial_sel[:, 3] == 4, trial_sel[:, 5] == 6)
        trial_sel_right_3 = np.logical_and(trial_sel[:, 3] == 8, trial_sel[:, 5] == 3)
        trial_sel_right_6 = np.logical_and(trial_sel[:, 3] == 8, trial_sel[:, 5] == 6)

        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(trial_FR[:, :, su_idx]).T

            left_trials_3 = onesu[trial_sel_left_3, :][:, self.row_sel_3]
            right_trials_3 = onesu[trial_sel_right_3, :][:, self.row_sel_3]

            left_trials_6 = onesu[trial_sel_left_6, :][:, self.row_sel_6]
            right_trials_6 = onesu[trial_sel_right_6, :][:, self.row_sel_6]

            left_trials = np.concatenate((left_trials_3, left_trials_6))

            right_trials = np.concatenate((right_trials_3, right_trials_6))

            bins = np.ones(56)

            for bin_idx in range(left_trials.shape[1]):
                try:
                    (stat, p) = stats.mannwhitneyu(
                        left_trials[:, bin_idx].flatten(),
                        right_trials[:, bin_idx].flatten(),
                        alternative="two-sided",
                    )
                    bins[bin_idx] = p < 0.05

                except ValueError:
                    bins[bin_idx] = 0

            self.statistical_test_selective.append(bins)

    def addSampleSelect(self, trial_FR, trial_sel):
        trial_sel_left = trial_sel[:, 2] == 4
        trial_sel_right = trial_sel[:, 2] == 8

        trial_sel_test_left = trial_sel[:, 3] == 4
        trial_sel_test_right = trial_sel[:, 3] == 8

        trial_sel_3 = trial_sel[:, 5] == 3
        trial_sel_6 = trial_sel[:, 5] == 6

        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(trial_FR[:, :, su_idx]).T
            (base_mean, base_std) = self.baselineVector(onesu)

            left_3_trials = onesu[np.logical_and(trial_sel_left, trial_sel_3), :][:, self.row_sel_3]
            left_6_trials = onesu[np.logical_and(trial_sel_left, trial_sel_6), :][:, self.row_sel_6]

            right_3_trials = onesu[np.logical_and(trial_sel_right, trial_sel_3), :][:, self.row_sel_3]
            right_6_trials = onesu[np.logical_and(trial_sel_right, trial_sel_6), :][:, self.row_sel_6]

            left_trials = np.concatenate((left_3_trials, left_6_trials))
            right_trials = np.concatenate((right_3_trials, right_6_trials))

            left_3_test_trials = onesu[np.logical_and(trial_sel_test_left, trial_sel_3), :][:, self.row_sel_3]
            left_6_test_trials = onesu[np.logical_and(trial_sel_test_left, trial_sel_6), :][:, self.row_sel_6]

            right_3_test_trials = onesu[np.logical_and(trial_sel_test_right, trial_sel_3), :][:, self.row_sel_3]
            right_6_test_trials = onesu[np.logical_and(trial_sel_test_right, trial_sel_6), :][:, self.row_sel_6]

            left_test_trials = np.concatenate((left_3_test_trials, left_6_test_trials))
            right_test_trials = np.concatenate((right_3_test_trials, right_6_test_trials))

            delay_3s_trials = onesu[trial_sel_3, :][:, self.row_sel_3]
            delay_6s_trials = onesu[trial_sel_6, :][:, self.row_sel_6]

            scaled_FR = np.concatenate((delay_3s_trials, delay_6s_trials))

            self.modulation.append(
                np.abs(np.mean((scaled_FR - base_mean) / base_std, axis=0))
            )

            self.heatmapL.append(np.mean((left_trials - base_mean) / base_std, axis=0))
            self.heatmapR.append(np.mean((right_trials - base_mean) / base_std, axis=0))

            self.test_heatmapL.append(np.mean((left_test_trials - base_mean) / base_std, axis=0))
            self.test_heatmapR.append(np.mean((right_test_trials - base_mean) / base_std, axis=0))

            self.selectivity.append(
                self.gauss_average(
                    (np.mean(right_trials, axis=0) - np.mean(left_trials, axis=0))
                    / (
                            np.mean(right_trials, axis=0)
                            + np.mean(left_trials, axis=0)
                            + 0.25
                    )
                )
            )

    def addPairSelect(self, trial_FR, trial_sel):

        trial_sel_left_3 = np.logical_and(
            trial_sel[:, 2] != trial_sel[:, 3], trial_sel[:, 5] == 3
        )
        trial_sel_right_3 = np.logical_and(
            trial_sel[:, 2] == trial_sel[:, 3], trial_sel[:, 5] == 3
        )

        trial_sel_left_6 = np.logical_and(
            trial_sel[:, 2] != trial_sel[:, 3], trial_sel[:, 5] == 6
        )
        trial_sel_right_6 = np.logical_and(
            trial_sel[:, 2] == trial_sel[:, 3], trial_sel[:, 5] == 6
        )

        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(trial_FR[:, :, su_idx]).T
            (base_mean, base_std) = self.baselineVector(onesu)
            left_trials_3 = onesu[
                            trial_sel_left_3, 20:52
                            ]  # 28 is test start, aka, -2s to +6s
            right_trials_3 = onesu[trial_sel_right_3, 20:52]

            left_trials_6 = onesu[trial_sel_left_6, 32:64]
            right_trials_6 = onesu[trial_sel_right_6, 32:64]

            left_trials = np.concatenate((left_trials_3, left_trials_6))
            right_trials = np.concatenate((right_trials_3, right_trials_6))

            self.pair_heatmapL.append(
                np.mean((left_trials - base_mean) / base_std, axis=0)
            )
            self.pair_heatmapR.append(
                np.mean((right_trials - base_mean) / base_std, axis=0)
            )

            self.pair_selectivity.append(
                self.gauss_average(
                    (np.mean(left_trials, axis=0) - np.mean(right_trials, axis=0))
                    / (
                            np.mean(right_trials, axis=0)
                            + np.mean(left_trials, axis=0)
                            + 0.25
                    )
                )
            )
            # pair stastic selective
            bins = np.ones(left_trials.shape[1])
            for bin_idx in range(left_trials.shape[1]):
                try:
                    (stat, p) = stats.mannwhitneyu(
                        left_trials[:, bin_idx].flatten(),
                        right_trials[:, bin_idx].flatten(),
                        alternative="two-sided",
                    )
                    bins[bin_idx] = p < 0.05
                except ValueError:
                    bins[bin_idx] = 0
            self.statistical_pair_selective.append(bins)

    def plotOneHeatmapPanel(self, data, ax, ylbl):
        im = plt.imshow(data, cmap="jet", aspect="auto", vmin=-2, vmax=2)
        [plt.plot([x, x], ax.get_ylim(), "-w") for x in np.array([1, 2]) * 4 - 0.5]

        if ylbl:
            ax.set_ylabel("Unit #")
        else:
            ax.set_yticklabels([])
        return im

    def plotSampleHeatmap(self, gs_outer_id, fh):
        ax_1 = fh.add_subplot(gs_outer_id)
        gs_1_L = gridspec.GridSpecFromSubplotSpec(10, 11, subplot_spec=gs_outer_id)

        L_arr = np.array(self.heatmapL)
        R_arr = np.array(self.heatmapR)

        sort_idx = np.argsort(
            np.mean(R_arr[:, 12:28], axis=1) - np.mean(L_arr[:, 12:28], axis=1)
        )

        ax_1_L = fh.add_subplot(gs_1_L[:, :5])
        self.plotOneHeatmapPanel(L_arr[sort_idx, 8:28], ax_1_L, True)
        ax_1_L.set_title("Sample1")

        ax_1_L.set_xticks(np.array([1, 5]) * 4 - 0.5)
        ax_1_L.set_xticklabels([0, 4])

        ax_1_R = fh.add_subplot(gs_1_L[:, 5:])
        im = self.plotOneHeatmapPanel(R_arr[sort_idx, 8:28], ax_1_R, False)
        ax_1_R.set_title("Sample2")
        ax_1_R.set_xticks(np.array([1, 5]) * 4 - 0.5)
        ax_1_R.set_xticklabels([0, 4])

        plt.colorbar(im, ticks=[-2, 0, 2], format="%d")

        ax_1.set_frame_on(False)
        ax_1.tick_params(axis="x", colors="w")
        ax_1.tick_params(axis="y", colors="w")
        ax_1.set_xlabel("Time (s)")

    def plotSampleSelectivity(self, gs_outer_id, fh):
        ax = fh.add_subplot(gs_outer_id)
        data = np.array(self.selectivity)[:, 8:28]
        sort_idx = np.flip(np.argsort(np.mean(data[:, 4:20], axis=1)))
        im = plt.imshow(
            data[sort_idx, :], cmap="jet", aspect="auto", vmin=-0.5, vmax=0.5
        )
        [plt.plot([x, x], ax.get_ylim(), "-w") for x in np.array([1, 2]) * 4 - 0.5]
        ax.set_xticks(np.array([1, 5]) * 4 - 0.5)
        ax.set_xticklabels([0, 4])
        ax.set_ylabel("Unit #")
        ax.set_title("Sample selectivity")
        ax.set_xlabel("Time (s)")
        plt.colorbar(im, ticks=[-0.5, 0, 0.5], format="%0.1f")

    def plotPairHeatmap(self, gs_outer_id, fh):
        ax_1 = fh.add_subplot(gs_outer_id)
        gs_1_L = gridspec.GridSpecFromSubplotSpec(10, 11, subplot_spec=gs_outer_id)

        L_arr = np.array(self.pair_heatmapL)
        R_arr = np.array(self.pair_heatmapR)

        sort_idx = np.argsort(
            np.mean(R_arr[:, 8:24], axis=1) - np.mean(L_arr[:, 8:24], axis=1)
        )

        ax_1_L = fh.add_subplot(gs_1_L[:, :5])
        self.plotOneHeatmapPanel(L_arr[sort_idx, 4:24], ax_1_L, True)
        ax_1_L.set_title("Pair")
        ax_1_L.set_xticks(np.array([1, 5]) * 4 - 0.5)
        ax_1_L.set_xticklabels(["T", "+4"])

        ax_1_R = fh.add_subplot(gs_1_L[:, 5:])
        im = self.plotOneHeatmapPanel(R_arr[sort_idx, 4:24], ax_1_R, False)
        ax_1_R.set_title("non-Pair")
        ax_1_R.set_xticks(np.array([1, 5]) * 4 - 0.5)
        ax_1_R.set_xticklabels(["T", "+4"])
        plt.colorbar(im, ticks=[-2, 0, 2], format="%d")
        ax_1.set_frame_on(False)
        ax_1.tick_params(axis="x", colors="w")
        ax_1.tick_params(axis="y", colors="w")
        ax_1.set_xlabel("Time (s)")

    def plotPairSelectivity(self, gs_outer_id, fh):
        ax = fh.add_subplot(gs_outer_id)
        data = np.array(self.pair_selectivity)
        sort_idx = np.flip(np.argsort(np.mean(data, axis=1)))
        im = plt.imshow(
            data[sort_idx, 4:28], cmap="jet", aspect="auto", vmin=-0.5, vmax=0.5
        )
        [plt.plot([x, x], ax.get_ylim(), "-w") for x in np.array([1, 2]) * 4 - 0.5]
        ax.set_xticks(np.array([1, 6]) * 4 - 0.5)
        ax.set_xticklabels(["T", "+5"])
        ax.set_ylabel("Unit #")
        ax.set_title("Pairing selectivity")
        ax.set_xlabel("Time (s)")
        plt.colorbar(im, ticks=[-0.5, 0, 0.5], format="%0.1f")

    def plotFracSampleSel(self, gs_outer_id, fh):
        ax = fh.add_subplot(gs_outer_id)
        ms = self.gauss_average(np.mean(self.statistical_sample_selective, axis=0))
        mt = self.gauss_average(np.mean(self.statistical_test_selective, axis=0))
        plt.plot(ms[8:52])
        plt.plot(mt[8:52])
        yspan = (0, 0.3)
        [
            plt.plot([x, x], yspan, "--", color="gray")
            for x in np.array([1, 2, 5, 6]) * 4 - 0.5
        ]
        ax.set_ylim(yspan)
        ax.set_xticks([3.5, 7.5, 19.5, 23.5])
        ax.set_xticklabels(["S", '+1', 'T', "+1"])
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Fraction of selective SU")
        plt.legend(['sample', 'test'])
        ax.set_title("Cue selectivity")

    def plotFracPairSel(self, gs_outer_id, fh):
        ax = fh.add_subplot(gs_outer_id)
        mm = self.gauss_average(np.mean(self.statistical_pair_selective, axis=0))

        plt.plot(mm[4:28])

        yspan = (0, 1)
        [
            plt.plot([x, x], yspan, "--", color="gray")
            for x in np.array([1, 2]) * 4 - 0.5
        ]
        ax.set_ylim(yspan)
        ax.set_xticks([3.5, 23.5])
        ax.set_xticklabels(["T", "+5"])
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Fraction of selective SU")
        ax.set_title("Pairing selectivity")

    def plotOverallModulation(self, gs_outer_id, fh):
        ax = fh.add_subplot(gs_outer_id)
        mm = self.gauss_average(np.mean(self.modulation, axis=0))

        plt.plot(mm[8:])

        yspan = (0, 3)
        [
            plt.plot([x, x], yspan, "--", color="gray")
            for x in np.array([1, 2, 5, 6]) * 4 - 0.5
        ]
        ax.set_ylim(yspan)
        ax.set_xticks([3.5, 7.5, 19.5, 23.5])
        ax.set_xticklabels(["S", "+1", "T", "+1"])
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Absolute normailzed firing rate")
        ax.set_title("Overall modulation")

    def plotFracHitMissSel(self, gs_outer_id, fh):
        ax = fh.add_subplot(gs_outer_id)
        mHM = self.gauss_average(np.mean(self.statistical_hitMiss_selective, axis=0))
        plt.plot(mHM[8:])
        if len(self.statistical_CRFalse_selective) > 0:
            mCF = self.gauss_average(np.mean(self.statistical_CRFalse_selective, axis=0))
            plt.plot(mCF[8:])
            plt.legend(['Hit-Miss', 'CR-False'])

        yspan = (0, 1)
        [
            plt.plot([x, x], yspan, "--", color="gray")
            for x in np.array([1, 2, 5, 6]) * 4 - 0.5
        ]
        ax.set_ylim(yspan)
        ax.set_xticks([3.5, 7.5, 19.5, 23.5])
        ax.set_xticklabels(["S", "+1", "T", "+1"])
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Fraction of selective SU")
        ax.set_title("Hit v.s. Miss")

    def plotAvgAuroc(self, gs_outer_id, fh):
        ax = fh.add_subplot(gs_outer_id)
        data = np.array(self.auroc)[:, 8:44]
        sort_idx = np.flip(np.argsort(np.mean(data[:, 4:20], axis=1)))
        im = plt.imshow(
            data[sort_idx, :], cmap="jet", aspect="auto", vmin=0.5, vmax=0.7
        )
        [plt.plot([x, x], ax.get_ylim(), "-w") for x in np.array([1, 2, 5, 6]) * 4 - 0.5]
        ax.set_xticks(np.array([1, 2, 5, 6]) * 4 - 0.5)
        ax.set_xticklabels(['S', '+1', 'T', '+1'])
        ax.set_ylabel("Unit #")
        ax.set_title("AUROC for sample")
        ax.set_xlabel("Time (s)")
        plt.colorbar(im, ticks=[0.5, 0.6, 0.7], format="%0.1f")

    def plotSummary(self):
        fh = plt.figure(figsize=[7.5, 10])
        gs_outer = gridspec.GridSpec(3, 3, figure=fh)

        self.plotSampleHeatmap(gs_outer[0], fh)
        self.plotSampleSelectivity(gs_outer[1], fh)

        self.plotFracSampleSel(gs_outer[2], fh)

        self.plotPairHeatmap(gs_outer[3], fh)
        self.plotPairSelectivity(gs_outer[4], fh)

        self.plotFracPairSel(gs_outer[5], fh)

        self.plotAvgAuroc(gs_outer[6], fh)

        self.plotOverallModulation(gs_outer[7], fh)

        self.plotFracHitMissSel(gs_outer[8], fh)

        fh.suptitle(self.regionName)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()
        fh.savefig(self.regionName + ".png", dpi=300, bbox_inches="tight")

    def plot_3_vs_6(self):
        # fh = plt.figure(figsize=[7.5, 10])
        # gs_outer = gridspec.GridSpec(3, 3, figure=fh)
        pass

    def getFeatureVector(self):
        return np.concatenate(
            (
                np.mean(self.statistical_sample_selective, axis=0)[8:44],
                np.mean(self.statistical_test_selective, axis=0)[24:44],
                np.mean(self.statistical_pair_selective, axis=0),
                np.mean(self.statistical_hitMiss_selective, axis=0)[8:],
                np.mean(self.statistical_CRFalse_selective, axis=0)[8:],
                np.mean(self.modulation[8:], axis=0)
            )
        )

    def getPerSUFeatureVector(self):

        return np.concatenate((self.heatmapL,
                               self.heatmapR,
                               self.test_heatmapL,
                               self.test_heatmapL,
                               self.pair_heatmapL,
                               self.pair_heatmapR,
                               self.modulation), axis=1)

    # 36, 20, 32 ,48, 48, 48
