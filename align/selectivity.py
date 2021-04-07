# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:37:14 2019

@author: Libra
"""
import csv
import os
import re
import h5py
import pandas as pd
import numpy as np
import align.su_region_align as align


# import scipy.stats as sp


def exact_mc_perm_test(x, xall, y, yall, nmc):
    k = 0
    diff = np.abs(x / xall - y / yall)
    for j in range(nmc):
        a = np.random.permutation(xall + yall)
        newx = np.sum(a[0: (x + y)] < xall)
        k += diff <= np.abs(newx / xall - (x + y - newx) / yall)
    return k / nmc


def plotOneBar(stats, label, filename):
    import matplotlib.pyplot as plt

    fh = plt.figure(1, figsize=(15, 2.5), dpi=300)
    XLim = 0
    maxy = 0
    bhMiss = []
    nCount = []
    statsStr = []
    for (idx, row) in enumerate(stats):
        if row[4] > 30:
            if row[3] > 5:
                ratio = row[1] / row[3]
                bhLearn = plt.bar(XLim + 1, ratio, width=1, color="r")
                plt.errorbar(XLim + 1, ratio, row[5][1] - ratio, fmt="-k.", lw=0.5)
                maxy = max(maxy, ratio)
            else:
                bhMiss = plt.bar(
                    XLim + 1,
                    1,
                    width=1,
                    color="w",
                    linestyle=":",
                    linewidth=0.5,
                    edgecolor="gray",
                )
            if row[4] > 5:
                ratio = row[2] / row[4]
                bhWT = plt.bar(XLim + 2, ratio, width=1, color="c")
                plt.errorbar(XLim + 2, ratio, row[6][1] - ratio, fmt="-k.", lw=0.5)
                maxy = max(maxy, ratio)
            else:
                bhMiss = plt.bar(
                    XLim + 2,
                    1,
                    width=1,
                    color="w",
                    linestyle=":",
                    linewidth=0.5,
                    edgecolor="gray",
                )

            if row[3] >= 100 and row[4] >= 100:

                p = exact_mc_perm_test(row[1], row[3], row[2], row[4], 1000)
                if p < 0.001:
                    statsStr.append([XLim + 1.5, "***"])
                elif p < 0.01:
                    statsStr.append([XLim + 1.5, "**"])
                elif p < 0.05:
                    statsStr.append([XLim + 1.5, "*"])
                else:
                    statsStr.append([XLim + 1.5, "NS"])
            # elif row[3] >= 50 and row[4] >= 50:
            #     statsStr.append([XLim + 1.5, ">50"])
            elif row[3] >= row[4]:
                statsStr.append([XLim + 1.5, "WT " + str(row[4])])
            else:
                statsStr.append([XLim + 1.5, "LN " + str(row[3])])
            #        plt.text(XLim+1,0.125,row[3],ha='center')
            #        plt.text(XLim+2,0.075,row[4],ha='center')
            nCount.append([idx] + row[3:5])
            XLim += 4
    for r in statsStr:
        if r[1].startswith("WT") or r[1].startswith("LN"):
            plt.text(
                r[0], maxy * 1.1, r[1], rotation=90, ha="center", va="top", color="gray"
            )
        else:
            plt.text(
                r[0], maxy * 1.1, r[1], rotation=90, ha="center", va="top", color="k"
            )
    plt.ylabel(label)
    ax = fh.gca()
    ax.set_xticks(
        np.linspace(1.5, (len(nCount) - 1) * 4 + 1.5, num=len(nCount), endpoint=True)
    )

    tics = []
    for row in nCount:
        tics.append(regionSet[row[0]])

    ax.set_xticklabels(tics, rotation=90)
    if not bhMiss:
        fh.legend(
            (bhLearn, bhWT, bhMiss),
            ("learning", "welltrained", "missing"),
            loc="upper right",
        )
    else:
        fh.legend((bhLearn, bhWT), ("learning", "welltrained"), loc="upper right")

    plt.ylim(0, maxy * 1.125)
    plt.xlim(0, len(nCount) * 4 + 1)
    plt.show()
    fh.savefig(filename, dpi=300, bbox_inches="tight")
    return nCount




def get_trials(path):
    trials = None
    if not os.path.isfile(os.path.join(path, "events.hdf5")):
        print("missing one events file, path is ", path)
        input("press Enter to continue")
    with h5py.File(os.path.join(path, "events.hdf5"), "r") as fe:
        dset = fe["trials"]
        trials = np.array(dset, dtype="int32")
    return trials


if __name__ == "__main__":
    # %% setup

    FR_Th = 1.0

    regionMatched = []
    regionL = align.getRegionList()
    # conversion = []

    # %% main loop
    for path in align.traverse("K:/neupix/DataSum/"):

        sampSel = []
        pairSel = []
        if not os.path.isfile(os.path.join(path, "selectivity.hdf5")):
            if not os.path.isfile(os.path.join(path, "NOSU")):
                print("missing one selectivity file, path is", path)
                input("press Enter to continue")
            continue

        with h5py.File(os.path.join(path, "selectivity.hdf5"), "r") as fe:
            dset = fe["sampleSel"]
            sampSel = np.transpose(np.array(dset, dtype="double"))
            dset = fe["pairSel"]
            pairSel = np.transpose(np.array(dset, dtype="double"))

        trials = get_trials(path)
        if trials is None:
            continue

        (perf_desc, perf_code, inWindow, correctResp) = align.judgePerformance(trials)


        (bs_id, time_s, who_did) = align.get_bsid_duration_who(path)
        (mice_id, date, imec_no) = get_miceid_date_imecno(path)

        if not (mice_id and date and imec_no):
            print("failed extraction metadata from ", path)
            input("press Enter to continue")
            continue

        depthL = align.getTrackRegion(regionL, mice_id, date, imec_no, who_did)

        unitInfo = pd.read_csv(os.path.join(path, "cluster_info.tsv"), sep="\t")

        spkNThresh = time_s * FR_Th

        good_su = unitInfo.loc[
            (unitInfo["KSLabel"] == "good") & (unitInfo["n_spikes"] >= spkNThresh),
            ["id", "depth"],
        ]

        for (index, row) in good_su.iterrows():
            suDepth = row["depth"]
            # fullR = matchDepth(suDepth, depthL)
            # reg = combineSubRegion(fullR)
            reg = align.matchDepth(suDepth, depthL, date, mice_id, imec_no)
            # conversion.append([fullR, reg])

            if np.isin(np.double(row["id"]), sampSel[:, 0]):
                regionMatched.append(
                    [path, row["id"], suDepth, reg, perf_code]
                    + sampSel[sampSel[:, 0] == np.double(row["id"]), 1:5].tolist()[0]
                    + pairSel[pairSel[:, 0] == np.double(row["id"]), 1:3].tolist()[0]
                    + [int(mice_id)]
                )

                if len(regionMatched) % 500 == 0:
                    print("processed ", len(regionMatched), " SUs")

            else:
                print("not in selectivity list :", path)

    allRegions = []
    allPerfType = []
    allPath = []
    for u in regionMatched:
        allRegions.append(u[3])
        allPerfType.append(u[4])
        allPath.append(u[0])

    regionSet = list(set(allRegions))
    regionSet.sort()
    pathSet = list(set(allPath))

    suMat = []
    for u in regionMatched:
        # pathIdx, regionIdx, id, depth,perfType
        suMat.append(
            [pathSet.index(u[0]), regionSet.index(u[3])]
            + u[1:2]
            + u[4:11]
            + u[2:3]
            + u[11:12]
        )
    suMat = np.array(suMat)
    np.save("SUSelMat.npy", suMat)

    sampStat = []
    delayStat = []
    testStat = []
    import statsmodels.stats.proportion as prop

    for idx in np.unique(suMat[:, 1]):
        learnStat = suMat[np.logical_and(suMat[:, 1] == idx, suMat[:, 3] == 2), 5]
        wtStat = suMat[np.logical_and(suMat[:, 1] == idx, suMat[:, 3] == 3), 5]
        sampStat.append(
            [
                regionSet[np.int32(idx)],
                np.sum(learnStat < 0.001),
                np.sum(wtStat < 0.001),
                learnStat.shape[0],
                wtStat.shape[0],
                prop.proportion_confint(
                    np.sum(learnStat < 0.001), learnStat.shape[0], method="normal"
                ),
                prop.proportion_confint(
                    np.sum(wtStat < 0.001), wtStat.shape[0], method="normal"
                ),
            ]
        )
        learnStat = suMat[np.logical_and(suMat[:, 1] == idx, suMat[:, 3] == 2), 7]
        wtStat = suMat[np.logical_and(suMat[:, 1] == idx, suMat[:, 3] == 3), 7]
        delayStat.append(
            [
                regionSet[np.int32(idx)],
                np.sum(learnStat < 0.001),
                np.sum(wtStat < 0.001),
                learnStat.shape[0],
                wtStat.shape[0],
                prop.proportion_confint(
                    np.sum(learnStat < 0.001), learnStat.shape[0], method="normal"
                ),
                prop.proportion_confint(
                    np.sum(wtStat < 0.001), wtStat.shape[0], method="normal"
                ),
            ]
        )
        learnStat = suMat[np.logical_and(suMat[:, 1] == idx, suMat[:, 3] == 2), 9]
        wtStat = suMat[np.logical_and(suMat[:, 1] == idx, suMat[:, 3] == 3), 9]
        testStat.append(
            [
                regionSet[np.int32(idx)],
                np.sum(learnStat < 0.001),
                np.sum(wtStat < 0.001),
                learnStat.shape[0],
                wtStat.shape[0],
                prop.proportion_confint(
                    np.sum(learnStat < 0.001), learnStat.shape[0], method="normal"
                ),
                prop.proportion_confint(
                    np.sum(wtStat < 0.001), wtStat.shape[0], method="normal"
                ),
            ]
        )

    plotOneBar(sampStat, "sample selective dur. sample", "sampleSel.png")
    plotOneBar(delayStat, "sample selective dur. delay", "delaySel.png")
    nCount = plotOneBar(testStat, "pair selective dur. test", "pairSel.png")
    for row in nCount:
        print(
            regionSet[row[0]], ", learning n = ", row[1], ", welltrained n = ", row[2]
        )

    # %% write csv for correlation with opto_gene
    with open("delayStats.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in delayStat:
            if row[3] >= 30 and row[4] >= 30:
                cwriter.writerow([row[0], row[1] / row[3], row[2] / row[4]])

    with open("testStats.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in testStat:
            if row[3] >= 30 and row[4] >= 30:
                cwriter.writerow([row[0], row[1] / row[3], row[2] / row[4]])

# count=np.zeros((len(regionSet),4))
#
# for idx in range(len(allRegions)):
#    count[regionSet.index(allRegions[idx]),allPerfType[idx]]+=1
#
# for r in range(len(regionSet)):
#    print('%s, %d, %d, %d' % (regionSet[r],count[r,0],count[r,2],count[r,3]))


#
#    os.chdir(path)
#    (goodU,mU)=countNeurons()
#    goodSum.append(goodU)
#    mUSum.append(mU)
#    paths.append(path)
#    return (goodSum,mUSum,paths)