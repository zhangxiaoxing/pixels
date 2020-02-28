# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:14:22 2020

@author: Libra
"""
import os
import h5py
import csv
import numpy as np
import scipy.stats as stats
import selectivity as zpy


class GLM_stats:
    def __init__(self):
        self.sample_sel_only_sample = None
        self.sample_sel_only_ED = None
        self.sample_sel_only_LD = None
        self.sample_sel_SP_ED = None
        self.sample_sel_ED_LD = None
        self.sample_sel_SP_ED_LD = None

        self.non_sel_mod_only_SP = None
        self.non_sel_mod_only_ED = None
        self.non_sel_mod_only_LD = None
        self.non_sel_mod_SP_ED = None
        self.non_sel_mod_ED_LD = None
        self.non_sel_mod_SP_ED_LD = None
        self.non_mod = None

        # self.row_sel_6 = np.concatenate(
        #     (np.arange(16), np.arange(16, 40, 2), np.arange(40, 68))
        # )
        # self.row_sel_3 = np.arange(56)

        self.use_ranksum = True

    def bool_stats_test(self, A, B):
        ### TODO alternative use of perm_test
        if self.use_ranksum:
            try:
                # flatten is used instead of mean to control for variance
                # ranksum is supposedly insensitive to length of data
                (stat, p) = stats.mannwhitneyu(
                    A.flatten(),
                    B.flatten(),
                    alternative="two-sided",
                )
                return p < 0.001

            except ValueError:
                return False

    def exact_mc_perm_test(self, xs, ys, nmc):
        n, k = len(xs), 0
        diff = np.abs(np.mean(xs) - np.mean(ys))
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
        return k / nmc

    def processGLMStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None):

        ### TODO: when variables are none

        # FR_scale = np.concatenate(
        #     (trial_FR[:, (trials[:, 5] == 3) & welltrain_window & correct_resp, :][self.row_sel_3, :, :],
        #      trial_FR[:, (trials[:, 5] == 6) & welltrain_window & correct_resp, :][self.row_sel_6, :, :]), axis=1)
        #
        # trial_perf_sel = np.concatenate((trials[(trials[:, 5] == 3) & welltrain_window & correct_resp, :],
        #                                  trials[(trials[:, 5] == 6) & welltrain_window & correct_resp, :]), axis=0)

        firing_rate_6 = trial_FR[:, (trials[:, 5] == 6) & welltrain_window & correct_resp, :]
        trial_perf_sel = trials[(trials[:, 5] == 6) & welltrain_window & correct_resp, :]

        trial_sel_left = trial_perf_sel[:, 2] == 4
        trial_sel_test_left = trial_perf_sel[:, 3] == 4
        pair_sel = (trial_perf_sel[:, 2] != trial_perf_sel[:, 3])

        self.sample_sel_only_sample = np.zeros((1, trial_FR.shape[2]))
        self.sample_sel_only_ED = np.zeros_like(self.sample_sel_only_sample)
        self.sample_sel_only_LD = np.zeros_like(self.sample_sel_only_sample)
        self.sample_sel_SP_ED = np.zeros_like(self.sample_sel_only_sample)
        self.sample_sel_ED_LD = np.zeros_like(self.sample_sel_only_sample)
        self.sample_sel_SP_ED_LD = np.zeros_like(self.sample_sel_only_sample)

        self.non_sel_mod_only_SP = np.zeros_like(self.sample_sel_only_sample)
        self.non_sel_mod_only_ED = np.zeros_like(self.sample_sel_only_sample)
        self.non_sel_mod_only_LD = np.zeros_like(self.sample_sel_only_sample)
        self.non_sel_mod_SP_ED = np.zeros_like(self.sample_sel_only_sample)
        self.non_sel_mod_ED_LD = np.zeros_like(self.sample_sel_only_sample)
        self.non_sel_mod_SP_ED_LD = np.zeros_like(self.sample_sel_only_sample)
        self.non_mod = np.zeros_like(self.sample_sel_only_sample)

        ### TODO: selective only during sample

        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(firing_rate_6[:, :, su_idx]).T
            left_trials = onesu[trial_sel_left, :]
            right_trials = onesu[~trial_sel_left, :]

            left_test_trials = onesu[trial_sel_test_left, :]
            right_test_trials = onesu[~trial_sel_test_left, :]

            # 3d delay trials sample 12:16, Delay 16:28, ED 18:22, LD 22:26, DM 28:36, reward 32:36
            # 6d delay trials sample 12:16, Delay 16:40, ED 18:22, LD 34:38, DM 40:48, reward 44:48

            BS = np.arange(0, 10)
            SP = np.arange(12, 16)
            ED = np.arange(18, 28)
            LD = np.arange(28, 38)

            sample_sel_SP = self.bool_stats_test(left_trials[:, SP], right_trials[:, SP])
            sample_sel_ED = self.bool_stats_test(left_trials[:, ED], right_trials[:, ED])
            sample_sel_LD = self.bool_stats_test(left_trials[:, LD], right_trials[:, LD])

            self.sample_sel_only_sample[0, su_idx] = sample_sel_SP and (not sample_sel_ED) and (not sample_sel_LD)
            self.sample_sel_only_ED[0, su_idx] = sample_sel_ED and (not sample_sel_SP) and (not sample_sel_LD)
            self.sample_sel_only_LD[0, su_idx] = sample_sel_LD and (not sample_sel_ED) and (not sample_sel_SP)
            self.sample_sel_SP_ED[0, su_idx] = sample_sel_SP and sample_sel_ED and (not sample_sel_LD)
            self.sample_sel_ED_LD[0, su_idx] = (not sample_sel_SP) and sample_sel_ED and sample_sel_LD
            self.sample_sel_SP_ED_LD[0, su_idx] = sample_sel_SP and sample_sel_ED and sample_sel_LD

            nsm_SP = ((not sample_sel_SP) and
                      self.bool_stats_test(onesu[:, BS], onesu[:, SP]))
            nsm_ED = ((not sample_sel_ED) and
                      self.bool_stats_test(onesu[:, BS], onesu[:, ED]))
            nsm_LD = ((not sample_sel_LD) and
                      self.bool_stats_test(onesu[:, BS], onesu[:, LD]))

            self.non_sel_mod_only_SP[0, su_idx] = nsm_SP and (not nsm_ED) and (not nsm_LD)
            self.non_sel_mod_only_ED[0, su_idx] = nsm_ED and (not nsm_LD) and (not nsm_SP)
            self.non_sel_mod_only_LD[0, su_idx] = nsm_LD and (not nsm_ED) and (not nsm_SP)

            self.non_sel_mod_SP_ED[0, su_idx] = nsm_ED and nsm_SP and (not nsm_LD)
            self.non_sel_mod_ED_LD[0, su_idx] = nsm_ED and nsm_LD and (not nsm_SP)
            self.non_sel_mod_SP_ED_LD[0, su_idx] = nsm_LD and nsm_ED and nsm_SP

            self.non_mod[0, su_idx] = not (self.bool_stats_test(onesu[:, BS], onesu[:, SP]) or
                                           self.bool_stats_test(onesu[:, BS], onesu[:, ED]) or
                                           self.bool_stats_test(onesu[:, BS], onesu[:, LD]))

            ### TODO: selective only during sample

    def getFeatures(self):
        return np.concatenate((self.sample_sel_only_sample,
                               self.sample_sel_only_ED,
                               self.sample_sel_only_LD,
                               self.sample_sel_SP_ED,
                               self.sample_sel_ED_LD,
                               self.sample_sel_SP_ED_LD,

                               self.non_sel_mod_only_SP,
                               self.non_sel_mod_only_ED,
                               self.non_sel_mod_only_LD,
                               self.non_sel_mod_SP_ED,
                               self.non_sel_mod_ED_LD,
                               self.non_sel_mod_SP_ED_LD,
                               self.non_mod))


### all brain region entry point
def prepare_GLM():
    curr_stats = GLM_stats()
    all_sess_list = []
    reg_list = []
    dpath = None
    if os.path.exists("/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"):
        dpath = "/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"
    else:
        dpath = r"K:\neupix\DataSum"
    for path in zpy.traverse(dpath):
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
        suid_reg = []
        with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
            l = list(csv.reader(csvfile))[1:]
            suid_reg = [list(i) for i in zip(*l)]

        (perf_desc, perf_code, welltrain_window, correct_resp) = zpy.judgePerformance(trials)

        if perf_code != 3:
            continue

        curr_stats.processGLMStats(trial_FR, trials, welltrain_window, correct_resp)

        onesession = curr_stats.getFeatures()
        all_sess_list.append(onesession)

        reg_list.extend(suid_reg[1])
    return (all_sess_list, reg_list)


def plot_features():
    # TODO: move figure plot here
    pass


# %% main
def process_all(denovo=False):
    all_sess_arr = None
    reg_arr = None
    if denovo:
        ### save raw data file
        (all_sess_list, reg_list) = prepare_GLM()
        all_sess_arr = np.concatenate(tuple(all_sess_list), axis=1)
        reg_arr = np.array(reg_list)
        np.savez_compressed('GLM_stats.npz', all_sess_arr=all_sess_arr, reg_arr=reg_arr)
    else:
        ### load back saved raw data
        if not os.path.isfile("GLM_stats.npz"):
            print('missing data file!')
            return
        np.load('GLM_stats.npz')
        fstr = np.load('GLM_stats.npz')
        reg_arr = fstr['reg_arr']
        all_sess_arr = fstr['all_sess_arr']

    reg_set = list(set(reg_arr.tolist()))

    plot_figure = False
    if plot_figure:
        reg_n = []
        for i in range(len(reg_set)):
            reg_n.append(f'{np.count_nonzero(reg_arr == reg_set[i])} SU {reg_set[i]}')
        feature_tag = ['sample selective during sample',
                       'sample selective during early delay',
                       'sample selective during late delay',
                       'sample selective ONLY during late delay',
                       'sample selective during decision making',
                       'test selective during decision making',
                       'pair selective during decision making',
                       'pair selective ONLY during reward',
                       'nonselective modulation during late delay',
                       'nonselective modulation during decision making',
                       'unmodulated']
        import matplotlib.pyplot as plt
        fh = plt.figure(figsize=(11.34, 40), dpi=300)


    su_factors = []
    su_factors.append(reg_set)

    reg_count = []
    for reg in reg_set:
        reg_count.append(np.count_nonzero(reg_arr == reg))
    su_factors.append(reg_count)

    for feature in range(all_sess_arr.shape[0]):
        per_region_su_factor = np.zeros_like(reg_set, dtype=np.float64)
        for i in range(len(reg_set)):
            per_region_su_factor[i] = np.mean(all_sess_arr[feature, reg_arr == reg_set[i]])
        su_factors.append(per_region_su_factor.tolist())

        if plot_figure:
            reg_idx = np.argsort(per_region_su_factor)
            ax = plt.subplot(11, 1, feature + 1)
            plt.bar(np.arange(len(reg_set)), per_region_su_factor[reg_idx])
            ax.set_xticks(np.arange(len(reg_set)))
            ax.set_xticklabels(np.array(reg_n)[reg_idx], rotation=90, ha='center', va='top', fontsize=7.5)
            ax.set_xlim(-1, len(reg_set))
            ax.set_ylabel(feature_tag[feature])

    if plot_figure:
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        fh.savefig("coding_feature.png", dpi=300, bbox_inches="tight")
        plt.show()

    su_factors = [list(i) for i in zip(*su_factors)]

    ### export csv for matlab GLM
    with open("glm_coding_features.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in su_factors:
            cwriter.writerow(row)

if __name__=="__main__":
    process_all(True)