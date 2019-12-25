# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:37:14 2019

@author: Libra
"""
import csv
import os
import re
import phylib.utils._misc as phyutil
import numpy as np
import h5py
import scipy.stats as sp

site_file="K:/neupix/meta/NP tracks coarsely labelled depth1219.csv"
unlabledRecord=[];
def traverse(path):
    for (basepath, dirs, files) in os.walk(path):
        if "cluster_info.tsv" in files:
            yield basepath


def imecNo2side(who, date, imecNo, mid):
    if date=='191130' and mid=='49':
        return 'L'
    
    if date=='191101' and mid=='26':
        return 'R'
    
    
    if who == "HRM" and (int(date)) >= 191028:
        if imecNo == "1":
            return "R"
        elif imecNo == "0":
            return "L"
    else:
        if imecNo == "1":
            return "L"
        elif imecNo == "0":
            return "R"

    print("Error parsing imec No")
    return "X"


def get_bsid_duration_who(path):
    bs_id=0
    time_s=0
    files = os.listdir(path)
    for f in files:
        if f.endswith("ap.meta"):
            with open(os.path.join(path, f), "r") as file:
                for line in file:
                    bs_id_grps = re.match("imDatBsc_sn=(\\d{1,3})", line)
                    if bs_id_grps:
                         bs_id=bs_id_grps.group(1)
                    fs_grps=re.match('fileSizeBytes=(\\d+)',line)
                    if fs_grps:
                        time_s=int(fs_grps.group(1))/385/2/30000
                        
    if bs_id == "350":
        who = "ZHA"
    elif bs_id == "142":
        who = "HRM"
    else:
        who='UNKNOWN'
        print("Unknown BS id!")
        input("press Enter to continue") 
                   
    return (bs_id,time_s,who)


def getTrackRegion(regionL, mid, date, imecNo, who):
    tMin = 0
    depthL = []
    for row in regionL:
        if (
            row[2] == '20'+date
            and row[0] == mid
            and row[3] == imecNo2side(who, date, imecNo, mid)

        ):
            depthL.append([row[6], int(row[11]), int(row[10])])
            tMin = min(tMin, int(row[11]))
            # tMin = min(tMin, int(row[10]))
    for r in depthL:
        r[1] -= tMin
        r[2] -= tMin

    return depthL


def matchDepth(depth, depthL):
    if len(depthL) > 0:
        for row in depthL:
            if row[1] <= depth < row[2]:
                return row[0]
    unlabledRecord.append([date, mice_id, imec_no, depth])
    return "Unlabeled"


def judgePerformance(trials):
    if trials.shape[0] >= 40:
        correctResp = np.bitwise_xor(trials[:, 4] == trials[:, 5], trials[:, 6] == 1)
        inWindow = np.zeros((trials.shape[0],), dtype="bool")
        i = 40
        while i < trials.shape[0]:
            #            if np.sum(correctResp[i-40:i])>=32:
            if np.sum(correctResp[i - 40 : i]) >= 32:
                inWindow[i - 40 : i] = 1
            i += 1
        if np.sum(inWindow) >= 40:  # Well Trained
            return ("wellTrained", 3, trials[inWindow, :])

        else:
            inWindow = np.zeros((trials.shape[0],), dtype="bool")
            licks = trials[:, 6] == 1
            i = 40
            while i < trials.shape[0]:
                if np.sum(licks[i - 40 : i]) >= 16:  # Learning
                    inWindow[i - 40 : i] = 1
                i += 1
            if np.sum(inWindow) >= 40:
                return ("learning", 2, trials[inWindow, :])
            elif np.sum(licks) <= trials.shape[0] // 10:  # Passive
                return ("passive", 0, trials)
            else:
                return ("transition", 1, trials)  # Unlabled


def combineSubRegion(r):
    if re.match("CA[13]", r):
        return r
    if re.match("([A-Za-z-]+)[1-6/]{0,3}[ab]{0,1}", r):
        g = re.match("([A-Za-z-]+)[1-6/]{0,3}[ab]{0,1}", r)
        return g.group(1)
    else:
        return r


def plotOneBar(stats, label, filename):
    import matplotlib.pyplot as plt

    fh = plt.figure(1, figsize=(15, 2.5), dpi=300)
    XLim = 0
    maxy = 0
    bhMiss = []
    nCount = []
    statsStr = []
    for row in stats:
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

            csq = sp.chisquare(
                [row[1], row[3] - row[1]], f_exp=[row[2], row[4] - row[2]]
            )
            if csq[1] < 0.001:
                statsStr.append([XLim + 1.5, "***"])
            elif csq[1] < 0.01:
                statsStr.append([XLim + 1.5, "**"])
            elif csq[1] < 0.05:
                statsStr.append([XLim + 1.5, "*"])
            else:
                statsStr.append([XLim + 1.5, "NS"])
        elif row[3]>=50 and row[4]>=50:
            statsStr.append([XLim+1.5, ">50"])
        else:
            statsStr.append([XLim + 1.5, "WIP"])

        #        plt.text(XLim+1,0.125,row[3],ha='center')
        #        plt.text(XLim+2,0.075,row[4],ha='center')
        nCount.append(row[3:5])
        XLim += 4
    for r in statsStr:
        if r[1] == "WIP" or r[1]==">50":
            plt.text(
                r[0],
                maxy * 1.05,
                r[1],
                rotation=90,
                ha="center",
                va="center",
                color="gray",
            )
        else:
            plt.text(
                r[0],
                maxy * 1.05,
                r[1],
                rotation=90,
                ha="center",
                va="center",
                color="k",
            )
    plt.ylabel(label)
    ax = fh.gca()
    ax.set_xticks(
        np.linspace(
            1.5, (len(regionSet) - 1) * 4 + 1.5, num=len(regionSet), endpoint=True
        )
    )
    ax.set_xticklabels(regionSet, rotation=90)
    if not bhMiss:
        fh.legend(
            (bhLearn, bhWT, bhMiss),
            ("learning", "welltrained", "missing"),
            loc="upper right",
        )
    else:
        fh.legend((bhLearn, bhWT), ("learning", "welltrained"), loc="upper right")

    plt.ylim(0, maxy * 1.125)
    plt.xlim(0, len(regionSet) * 4 + 1)
    plt.show()
    fh.savefig(filename, dpi=300, bbox_inches="tight")
    return nCount

def get_region_list():
    regionL = []
    
    with open(site_file, newline="") as cf:
        creader = csv.reader(cf, dialect="excel")
        for row in creader:
            regionL.append(row)
    return regionL


def get_miceid_date_imecno(path):
    dateA = re.compile("(19\\d\\d\\d\\d)\\D")
    imecNo = re.compile("imec(\\d)")
    miceId = re.compile("[M\\\\\\-_](\\d{2})[\\\\_\\-]")

    dateGrps = dateA.search(path)
    imecGrps = imecNo.search(path)
    miceGrps=miceId.search(path)
    
    if dateGrps and imecGrps and miceGrps:
        # print(miceGrps.group(1),', ',path)
        return (miceGrps.group(1),dateGrps.group(1),imecGrps.group(1))
    else:
        print('unresolved path meta data: ',path)
        input('press Enter to continue')
        return (None,None,None)
    
def get_trials(path):
    trials = None
    if not os.path.isfile(os.path.join(path, "events.hdf5")):
        print("missing one events file, path is ", path)
        input("press Enter to continue")
    with h5py.File(os.path.join(path, "events.hdf5"), "r") as fe:
        dset = fe["trials"]
        trials = np.array(dset, dtype="int32")    
    return trials


if __name__=="__main__":
#%% setup
    

    
    FR_Th = 1.0
    

    
    regionMatched = []
    conversion = []
    
    # %% main loop
    for path in traverse("K:/neupix/DataSum/"):
    

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
    
        trials=get_trials(path)
        if trials is None:
            continue
        
        (perfType, perfIdx, selTrials) = judgePerformance(trials)
        
        (bs_id,time_s,who)=get_bsid_duration_who(path)
        (mice_id, date, imec_no)=get_miceid_date_imecno(path)
        
        regionL=get_region_list();
    
        depthL=getTrackRegion(regionL, mice_id, date, imec_no, who)
        
        unitInfo = phyutil.read_tsv(os.path.join(path, "cluster_info.tsv"))
        spkNThresh = time_s * FR_Th
    
        for row in unitInfo:
            if row["KSLabel"] == "good" and row["n_spikes"] >= spkNThresh:
                suDepth = row["depth"]
                fullR = matchDepth(suDepth, depthL)
                reg = combineSubRegion(fullR)
                conversion.append([fullR, reg])
    
                if np.isin(np.double(row["id"]), sampSel[:, 0]):
                    regionMatched.append(
                        [path, row["id"], suDepth, reg, perfIdx]
                        + sampSel[sampSel[:, 0] == np.double(row["id"]), 1:5].tolist()[0]
                        + pairSel[pairSel[:, 0] == np.double(row["id"]), 1:3].tolist()[0]
                        + [int(mice_id)]
                    )
                else:
                    print("not in selectivity list")
                    
    
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
    for idx in range(len(regionSet)):
        print(
            regionSet[idx],
            ", learning n = ",
            nCount[idx][0],
            ", welltrained n = ",
            nCount[idx][1],
        )



#%% write csv for correlation with opto_gene
with open('delayStats.csv','w',newline='') as cf:
    cwriter = csv.writer(cf, dialect="excel")
    for row in delayStat:
        if row[3]>=10 and row[4]>=10:
            cwriter.writerow([row[0],row[1]/row[3],row[2]/row[4]])        





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
