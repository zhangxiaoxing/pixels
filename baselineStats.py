# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:14:22 2020

@author: Libra
"""
import os
import h5py
import csv
import numpy as np
# from scipy.stats.stats import pearsonr
import scipy.stats as stats
import su_region_align as align


class baseline_stats:
    def __init__(self):
        # self.per_trial_base_avg = None
        self.stats = None
        # self.correct_base_avg = None
        # self.err_base_avg = None

    def process_baseline_stats(self, trial_FR, trials, welltrain_window=None, correct_resp=None):
        ### TODO: when variables are none

        # FR_scale = np.concatenate(
        #     (trial_FR[:, (trials[:, 5] == 3) & welltrain_window & correct_resp, :][self.row_sel_3, :, :],
        #      trial_FR[:, (trials[:, 5] == 6) & welltrain_window & correct_resp, :][self.row_sel_6, :, :]), axis=1)
        #
        # trial_perf_sel = np.concatenate((trials[(trials[:, 5] == 3) & welltrain_window & correct_resp, :],
        #                                  trials[(trials[:, 5] == 6) & welltrain_window & correct_resp, :]), axis=0)

        all_baseline = trial_FR[:12, :, :]
        correct_baseline = trial_FR[:, welltrain_window & correct_resp, :][:12, :, :]
        error_baseline = trial_FR[:, np.logical_not(correct_resp), :][:12, :, :]

        ### TODO: selective only during sample
        self.stats = []
        for su_idx in range(trial_FR.shape[2]):
            # print(su_idx)
            avg_fr = np.mean(all_baseline[:, :, su_idx], axis=0)
            onesu = np.vstack((np.arange(avg_fr.shape[0]), avg_fr)).T
            if np.unique(onesu[:, 1]).size < 3:
                (r, pr) = (0, 1)
            else:
                (r, pr) = stats.pearsonr(onesu[:, 0], onesu[:, 1])
            cr = np.mean(correct_baseline[:, :, su_idx], axis=0)
            err = np.mean(error_baseline[:, :, su_idx], axis=0)
            if cr.size > 30 and err.size > 30 and np.unique(np.hstack((cr, err))).size >= 2:
                p = stats.mannwhitneyu(cr, err, alternative='two-sided')[1]
            else:
                p = 1

            # self.per_trial_base_avg.append(onesu)
            self.stats.append([r, pr, p])

            ### TODO: selective only during sample

    def getFeatures(self):
        return self.stats


### all brain region entry point
def prepare_baseline():
    curr_stats = baseline_stats()
    all_sess_list = []
    reg_list = []
    paths = []
    suids = []
    dpath = align.get_root_path()
    # counter = 0
    for path in align.traverse(dpath):
        # path = r'K:\neupix\DataSum\200104_75_learning4_g1\200104_75_learning4_g1_imec0_cleaned'
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
                    dset = ffr["SU_id"]
                    SU_ids = np.array(dset, dtype="uint16")
                    dset = ffr["FR_All"]
                    trial_FR = np.array(dset, dtype="double")  # bin/68, trials/230, su/NNN
                    dset = ffr["Trials"]
                    trials = np.array(dset, dtype="double").T
                done_read = True
            except OSError:
                print("h5py read error handled")
        if trials is None:
            continue
        suid_reg = None
        with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
            l = list(csv.reader(csvfile))[1:]
            suid_reg = [list(i) for i in zip(*l)]

        (perf_desc, perf_code, welltrain_window, correct_resp) = align.judgePerformance(trials)

        if perf_code != 3:
            continue

        curr_stats.process_baseline_stats(trial_FR, trials, welltrain_window, correct_resp)

        onesession = curr_stats.getFeatures()
        all_sess_list.extend(onesession)
        reg_list.extend(suid_reg[1])
        for one_su in SU_ids.flatten():
            paths.append(path)
            suids.append(one_su)
        # counter += 1
        # if counter > 2:
        #     break

    export = []
    export.append(paths)
    export.append(suids)
    export.append(reg_list)
    export.extend(list(zip(*all_sess_list)))

    export = list(zip(*export))
    with open("baseline_stats.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in export:
            cwriter.writerow(row)

    return (all_sess_list, reg_list, paths, suids)


def plot_features():
    # TODO: move figure plot here
    pass


# %% main
def process_all(denovo=False):
    all_sess_arr = None
    reg_arr = None
    if denovo:
        ### save raw data file
        (all_sess_list, reg_list, paths, suids) = prepare_baseline()
        all_sess_arr = np.array(all_sess_list)
        reg_arr = np.array(reg_list)
        np.savez_compressed('baseline_stats.npz', all_sess_arr=all_sess_arr, reg_arr=reg_arr, paths=paths, suids=suids)
    else:
        ### load back saved raw data
        if not os.path.isfile("baseline_stats.npz"):
            print('missing data file!')
            return
        fstr = np.load('baseline_stats.npz')
        reg_arr = fstr['reg_arr']
        all_sess_arr = fstr['all_sess_arr']

    reg_set = list(set(reg_arr.tolist()))


if __name__ == "__main__":
    process_all(True)
