# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:49:17 2020

@author: Libra
"""
import os

import sys
import csv
import h5py
import numpy as np
# import scipy.stats as stats
# from sklearn.svm import SVC
from sklearn.model_selection import KFold
# from mcc.mcc import MaximumCorrelationClassifier
# from sklearn.model_selection import cross_val_score
# from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import su_region_align as align
# import per_second_stats
# from multiprocessing import Pool
# from matplotlib import rcParams
from sklearn.decomposition import PCA
from mpl_toolkits import mplot3d
import scipy.stats as stats
from pixelStats import baselineVector


class GLM_PCA_stats:
    def __init__(self):

        self.S1T1 = None
        self.S1T2 = None
        self.S2T1 = None
        self.S2T2 = None

        self.row_sel_6 = np.concatenate(
            (np.arange(16), np.arange(16, 40, 2), np.arange(40, 68))
        )
        self.row_sel_3 = np.arange(56)

    def processGLMStats(self, trial_FR, trials, welltrain_window=[], correct_resp=[]):

        ### TODO: when variables are empty

        FR_scale = np.concatenate(
            (trial_FR[:, (trials[:, 5] == 3) & welltrain_window & correct_resp, :][self.row_sel_3, :, :],
             trial_FR[:, (trials[:, 5] == 6) & welltrain_window & correct_resp, :][self.row_sel_6, :, :]), axis=1)

        trial_perf_sel = np.concatenate((trials[(trials[:, 5] == 3) & welltrain_window & correct_resp, :],
                                         trials[(trials[:, 5] == 6) & welltrain_window & correct_resp, :]), axis=0)

        # S1=4 S2=8 T1=8 T2 =4
        trial_S1T1 = (trial_perf_sel[:, 2] == 4) & (trial_perf_sel[:, 3] == 8)
        trial_S1T2 = (trial_perf_sel[:, 2] == 4) & (trial_perf_sel[:, 3] == 4)
        trial_S2T1 = (trial_perf_sel[:, 2] == 8) & (trial_perf_sel[:, 3] == 8)
        trial_S2T2 = (trial_perf_sel[:, 2] == 8) & (trial_perf_sel[:, 3] == 4)

        self.S1T1 = np.zeros((56, trial_FR.shape[2]))
        self.S1T2 = np.zeros((56, trial_FR.shape[2]))
        self.S2T1 = np.zeros((56, trial_FR.shape[2]))
        self.S2T2 = np.zeros((56, trial_FR.shape[2]))

        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(FR_scale[:, :, su_idx]).T

            (base_mean, base_std) = baselineVector(onesu)

            onesu = (onesu - base_mean) / base_std
            self.S1T1[:, su_idx] = np.mean(onesu[trial_S1T1, :], axis=0)
            self.S1T2[:, su_idx] = np.mean(onesu[trial_S1T2, :], axis=0)
            self.S2T1[:, su_idx] = np.mean(onesu[trial_S2T1, :], axis=0)
            self.S2T2[:, su_idx] = np.mean(onesu[trial_S2T2, :], axis=0)

    def getFeatures(self):

        return np.concatenate((self.S1T1,
                               self.S1T2,
                               self.S2T1,
                               self.S2T2))


class ctd_stats:
    def __init__(self):
        self.su_list = []
        self.kf = KFold(2, shuffle=True)
        self.avail = []


    def splitMean(self, FR, su_idx, mm, std, repeats=20):
        h1 = []
        h2 = []
        for i in range(repeats):
            sel = list(self.kf.split(np.arange(FR.shape[1])))[0]
            h1.append(np.squeeze(np.mean((FR[:, sel[0], su_idx] - mm) / std, axis=1)))
            h2.append(np.squeeze(np.mean((FR[:, sel[1], su_idx] - mm) / std, axis=1)))
        return (h1, h2)

        ### unprocessed data entry point

    def processCTDStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None):
        shifted_trial = np.vstack((np.zeros(6), trials[:-1, :]))

        ### TODO: when variables are empty

        FR_D3_S1_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 4, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D3_S2_S1 = trial_FR[
                      :, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 8, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :,
                      ]
        FR_D3_S1_S2 = trial_FR[
                      :, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 4, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :,
                      ]
        FR_D3_S2_S2 = trial_FR[
                      :, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 8, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :,
                      ]

        FR_D6_S1_S1 = trial_FR[
                      :, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 4, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :,
                      ]
        FR_D6_S2_S1 = trial_FR[
                      :, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 8, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :,
                      ]

        FR_D6_S1_S2 = trial_FR[
                      :, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 4, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :,
                      ]
        FR_D6_S2_S2 = trial_FR[
                      :, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 8, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :,
                      ]

        for su_idx in range(trial_FR.shape[2]):
            (mm, std) = baselineVector(
                np.squeeze(trial_FR[:, np.logical_and(welltrain_window, correct_resp), su_idx]))

            if std > 0 and np.all([x.shape[1] >= 2 for x in (
                    FR_D3_S1_S1,
                    FR_D3_S1_S2,
                    FR_D3_S2_S1,
                    FR_D3_S2_S2,
                    FR_D6_S1_S1,
                    FR_D6_S1_S2,
                    FR_D6_S2_S1,
                    FR_D6_S2_S2,
            )]):
                self.su_list.append(
                    {

                        "S1_S1_3m": self.splitMean(FR_D3_S1_S1, su_idx, mm, std),
                        "S1_S2_3m": self.splitMean(FR_D3_S1_S2, su_idx, mm, std),
                        "S2_S1_3m": self.splitMean(FR_D3_S2_S1, su_idx, mm, std),
                        "S2_S2_3m": self.splitMean(FR_D3_S2_S2, su_idx, mm, std),
                        "S1_S1_6m": self.splitMean(FR_D6_S1_S1, su_idx, mm, std),
                        "S1_S2_6m": self.splitMean(FR_D6_S1_S2, su_idx, mm, std),
                        "S2_S1_6m": self.splitMean(FR_D6_S2_S1, su_idx, mm, std),
                        "S2_S2_6m": self.splitMean(FR_D6_S2_S2, su_idx, mm, std),

                    }
                )
                self.avail.append(True)
            else:
                self.avail.append(False)

    def get_features(self):
        # breakpoint()
        return (self.su_list, self.avail)


def get_dataset(denovo):
    features_per_su = []
    avails = []
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
            (feat, avail) = currStats.get_features()
            features_per_su.extend(feat)
            avails.extend(avail)
            reg_list.extend(suid_reg[1])

            ### DEBUG in small subsets
            # if len(features_per_su)>50:
            #     break

        ### save to npz file
        np.savez_compressed("ctd_prev.npz", features_per_su=features_per_su, reg_list=reg_list, avails=avails)

    ### load from npz file
    else:
        fstr = np.load("ctd_prev.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        reg_list = fstr["reg_list"].tolist()
        avails = fstr['avaisl']
    return (features_per_su, reg_list, avails)


def plotPCA3D_FR(fr):
    pcamat = np.hstack([np.vstack([x['S1_S1_3m'] for x in fr]),
                        np.vstack([x['S2_S1_3m'] for x in fr]),
                        np.vstack([x['S1_S2_3m'] for x in fr]),
                        np.vstack([x['S2_S2_3m'] for x in fr]),
                        np.vstack([x['S1_S1_6m'] for x in fr]),
                        np.vstack([x['S2_S1_6m'] for x in fr]),
                        np.vstack([x['S1_S2_6m'] for x in fr]),
                        np.vstack([x['S2_S2_6m'] for x in fr]), ])
    np.save(f'pcamat_{delay}_sust.npy', pcamat)

    pca = PCA(n_components=20)
    comp = pca.fit_transform(pcamat.T)
    coeff = pca.components_
    ratio = pca.explained_variance_ratio_

    fig = plt.figure(figsize=(5, 5), dpi=300)
    ax = plt.axes(projection="3d")

    ax.plot3D(comp[(68 * 4):(68 * 5), 0], comp[(68 * 4):(68 * 5), 1], comp[(68 * 4):(68 * 5), 2], 'r-')  # S1S1
    ax.plot3D(comp[(68 * 5):(68 * 6), 0], comp[(68 * 5):(68 * 6), 1], comp[(68 * 5):(68 * 6), 2], 'b-')  # S2S1
    ax.plot3D(comp[(68 * 6):(68 * 7), 0], comp[(68 * 6):(68 * 7), 1], comp[(68 * 6):(68 * 7), 2], 'm-')
    ax.plot3D(comp[(68 * 7):(68 * 8), 0], comp[(68 * 7):(68 * 8), 1], comp[(68 * 7):(68 * 8), 2], 'c-')
    plt.show()
    return pcamat


if __name__ == "__main__":
    delay = 6
    (features_per_su, reg_list, avail) = get_dataset(True)
    trans_fstr = np.load(f'sus_trans_pie_{delay}.npz')
    # list(trans6.keys())
    sust = trans_fstr['sust']
    trans = trans_fstr['transient']
    # reg_arr = trans_fstr['reg_arr']
    
    trans_avail=trans[avail]
    
    fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
    # fr_sust = [f for (f, t) in zip(features_per_su, sust) if t]
    # plotPCA3D_FR(fr_sust)
    # plotPCA3D_FR(fr_trans)
    
    fr=fr_trans
    same_prev=np.zeros((20,68))
    oppo_prev=np.zeros((20,68))
    
    for bidx in range(fr[0]['S1_S1_6m'][0][0].size): # 68 bins
        for ridx in range(len(fr[0]['S1_S1_6m'][0])):
            S1S1H1=np.array([x['S1_S1_6m'][0][ridx][bidx] for x in fr])
            S1S1H2=np.array([x['S1_S1_6m'][1][ridx][bidx] for x in fr])
            S2S1H1=np.array([x['S2_S1_6m'][0][ridx][bidx] for x in fr])
            S2S1H2=np.array([x['S2_S1_6m'][1][ridx][bidx] for x in fr])

            S1S2H1=np.array([x['S1_S2_6m'][0][ridx][bidx] for x in fr])
            S1S2H2=np.array([x['S1_S2_6m'][1][ridx][bidx] for x in fr])
            S2S2H1=np.array([x['S2_S2_6m'][0][ridx][bidx] for x in fr])
            S2S2H2=np.array([x['S2_S2_6m'][1][ridx][bidx] for x in fr])


            
            S1S1In=np.sqrt(np.sum(np.square(S1S1H1-S1S1H2)))
            S1S2In=np.sqrt(np.sum(np.square(S1S2H1-S1S2H2)))
            S2S1In=np.sqrt(np.sum(np.square(S2S1H1-S2S1H2)))
            S2S2In=np.sqrt(np.sum(np.square(S2S2H1-S2S2H2)))
            
            S1Out1=np.sqrt(np.sum(np.square(S1S1H1-S2S1H1)))
            S1Out2=np.sqrt(np.sum(np.square(S1S1H1-S2S1H2)))
            S1Out3=np.sqrt(np.sum(np.square(S1S1H2-S2S1H1)))
            S1Out4=np.sqrt(np.sum(np.square(S1S1H2-S2S1H2)))
            
            S2Out1=np.sqrt(np.sum(np.square(S1S2H1-S2S2H1)))
            S2Out2=np.sqrt(np.sum(np.square(S1S2H1-S2S2H2)))
            S2Out3=np.sqrt(np.sum(np.square(S1S2H2-S2S2H1)))
            S2Out4=np.sqrt(np.sum(np.square(S1S2H2-S2S2H2)))
            
            same_prev[ridx,bidx]=np.mean([S1S1In,S1S2In,S2S1In,S2S2In])
            oppo_prev[ridx,bidx]=np.mean([S1Out1,S1Out2,S1Out3,S1Out4,S2Out1,S2Out2,S2Out3,S2Out4,])



    same_mm=np.mean(same_prev,axis=0)            
    oppo_mm=np.mean(oppo_prev,axis=0)
    same_sem=stats.sem(same_prev)
    oppo_sem=stats.sem(oppo_prev)
    
    wrs=np.zeros(same_mm.shape[0])
    for bidx in range(same_mm.shape[0]):
        wrs[bidx]=stats.ranksums(same_prev[:,bidx],oppo_prev[:,bidx])[1]

    (fig,ax)=plt.subplots(1,1,figsize=(2,2),dpi=300)
    plt.fill_between(np.arange(same_mm.shape[0]), same_mm - same_sem, same_mm + same_sem, color="b", alpha=0.2)
    plt.fill_between(np.arange(oppo_mm.shape[0]), oppo_mm - oppo_sem, oppo_mm + oppo_sem, color="r", alpha=0.2)
    
    ph0=ax.plot(same_mm,'b-',lw=1)
    ph1=ax.plot(oppo_mm,'r-',lw=1)
    for bidx in range(same_mm.shape[0]):
        if wrs[bidx]*same_mm.size < 0.05:
            ax.plot(bidx,120,'k.',markersize=1)
    
    ax.set_xticks([12,32,52])
    ax.set_xticklabels([0,5,10])
    ax.set_xlabel('time (s)')
    ax.set_ylabel('trajectory distance (A.U.)')
    [ax.axvline(x, color='k', ls=':',lw=0.5) for x in [11.5, 15.5, 39.5, 43.5]]
    ax.set_xlim([4,40])
    # ax.legend((ph0[0], ph1[0]), ('same','oppo.'))
    fig.savefig('transient_6_hist_traj_dist.pdf',bbox_inches='tight')

    plt.show()
    

    sust_avail=sust[avail]
    fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]

    fr=fr_sust
    same_prev=np.zeros((20,68))
    oppo_prev=np.zeros((20,68))
    
    for bidx in range(fr[0]['S1_S1_6m'][0][0].size): # 68 bins
        for ridx in range(len(fr[0]['S1_S1_6m'][0])):
            S1S1H1=np.array([x['S1_S1_6m'][0][ridx][bidx] for x in fr])
            S1S1H2=np.array([x['S1_S1_6m'][1][ridx][bidx] for x in fr])
            S2S1H1=np.array([x['S2_S1_6m'][0][ridx][bidx] for x in fr])
            S2S1H2=np.array([x['S2_S1_6m'][1][ridx][bidx] for x in fr])

            S1S2H1=np.array([x['S1_S2_6m'][0][ridx][bidx] for x in fr])
            S1S2H2=np.array([x['S1_S2_6m'][1][ridx][bidx] for x in fr])
            S2S2H1=np.array([x['S2_S2_6m'][0][ridx][bidx] for x in fr])
            S2S2H2=np.array([x['S2_S2_6m'][1][ridx][bidx] for x in fr])


            
            S1S1In=np.sqrt(np.sum(np.square(S1S1H1-S1S1H2)))
            S1S2In=np.sqrt(np.sum(np.square(S1S2H1-S1S2H2)))
            S2S1In=np.sqrt(np.sum(np.square(S2S1H1-S2S1H2)))
            S2S2In=np.sqrt(np.sum(np.square(S2S2H1-S2S2H2)))
            
            S1Out1=np.sqrt(np.sum(np.square(S1S1H1-S2S1H1)))
            S1Out2=np.sqrt(np.sum(np.square(S1S1H1-S2S1H2)))
            S1Out3=np.sqrt(np.sum(np.square(S1S1H2-S2S1H1)))
            S1Out4=np.sqrt(np.sum(np.square(S1S1H2-S2S1H2)))
            
            S2Out1=np.sqrt(np.sum(np.square(S1S2H1-S2S2H1)))
            S2Out2=np.sqrt(np.sum(np.square(S1S2H1-S2S2H2)))
            S2Out3=np.sqrt(np.sum(np.square(S1S2H2-S2S2H1)))
            S2Out4=np.sqrt(np.sum(np.square(S1S2H2-S2S2H2)))
            
            same_prev[ridx,bidx]=np.mean([S1S1In,S1S2In,S2S1In,S2S2In])
            oppo_prev[ridx,bidx]=np.mean([S1Out1,S1Out2,S1Out3,S1Out4,S2Out1,S2Out2,S2Out3,S2Out4,])



    same_mm=np.mean(same_prev,axis=0)            
    oppo_mm=np.mean(oppo_prev,axis=0)
    same_sem=stats.sem(same_prev)
    oppo_sem=stats.sem(oppo_prev)
    
    wrs=np.zeros(same_mm.shape[0])
    for bidx in range(same_mm.shape[0]):
        wrs[bidx]=stats.ranksums(same_prev[:,bidx],oppo_prev[:,bidx])[1]

    (fig,ax)=plt.subplots(1,1,figsize=(2,2),dpi=300)
    plt.fill_between(np.arange(same_mm.shape[0]), same_mm - same_sem, same_mm + same_sem, color="b", alpha=0.2)
    plt.fill_between(np.arange(oppo_mm.shape[0]), oppo_mm - oppo_sem, oppo_mm + oppo_sem, color="r", alpha=0.2)
    
    ph0=ax.plot(same_mm,'b-',lw=1)
    ph1=ax.plot(oppo_mm,'r-',lw=1)
    for bidx in range(same_mm.shape[0]):
        if wrs[bidx]*same_mm.size < 0.05:
            ax.plot(bidx,45,'k.',markersize=1)
    
    ax.set_xticks([12,32,52])
    ax.set_xticklabels([0,5,10])
    ax.set_xlabel('time (s)')
    ax.set_ylabel('trajectory distance (A.U.)')
    [ax.axvline(x, color='k', ls=':',lw=0.5) for x in [11.5, 15.5, 39.5, 43.5]]
    # ax.legend((ph0[0], ph1[0]), ('same','oppo.'))
    ax.set_xlim([4,40])
    fig.savefig('sust_6_hist_traj_dist.pdf',bbox_inches='tight')

    plt.show()
        
    
    
    
    
    
    
    
    