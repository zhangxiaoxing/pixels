# -*- coding: utf-8 -*-
"""
@author: Libra
"""

# import glob
import sys, os, glob
import numpy as np
from collections import deque
from readsync import readsync


def parseGNGEvents(events):
    """
    for Go/Nogo task event parsing
    """
    s1s = 30000  # sample per seconds
    trials = []
    lastTS = -300000
    lastCue = -1
    cueTS = -1
    rsps = -1
    cue = -1

    for eidx in range(len(events)):
        cue = events[eidx][1] & 0x0C  # bit mask S1, S2
        cueTS = events[eidx][0]  # Time stamp
        if cue > 0 and cueTS > (lastTS + s1s):  # update trial
            if lastCue > 0:
                trials.append([lastTS, lastCue, rsps])
                rsps = -1
            lastCue = cue
            lastTS = cueTS

        if (
            lastCue > 0
            and rsps < 0
            and events[eidx][0] >= lastTS + 30000
            and events[eidx][0] < lastTS + 90000
            and (events[eidx][1] & 0x03) > 0  # bit mask lick L/R
        ):
            rsps = events[eidx][1] & 0x03

    return trials


# def parseDPAEvents(events):
#     """
#     for DPA task parsing
#     if no cue,  assume sample
#     if has cue history check time diff
#     time diff < 3 discard
#     time diff > 3 < 10, assume current cue as test, previous cue as sample
#     time diff >10
# 
#     if sample and test assigned
#         if responsed
#             add responsed trial
#         else
#             add unresponsed trial
# 
#         reinit sample test response
#         assume current cue as sample"""
#     s1s = 30000
#     trials = []
#     lastTS = -300000 * 20
#     lastCue = -1
#     cueTS = -1
#     rsps = -1
#     cue = -1
#     sample = -1
#     test = -1
#     sampleTS = -1
#     testTS = -1
# 
#     for eidx in range(len(events)):
#         cue = events[eidx][1] & 0x0C  # bit mask S1, S2
#         cueTS = events[eidx][0]  # Time stamp
#         if cue > 0 and cueTS > lastTS + 1.1 * s1s and cueTS < lastTS + 2 * s1s:
#             print("error processing evt idx ", eidx)
#         elif (
#             cue > 0 and cueTS > lastTS + s1s * 2 and cueTS < lastTS + s1s * 8
#         ):  # following delay
#             sample = lastCue
#             sampleTS = lastTS
#             test = cue
#             testTS = cueTS
# 
#             lastCue = cue
#             lastTS = cueTS
# 
#         elif cue > 0 and cueTS > lastTS + s1s * 8:  # following ITI
#             if sample > 0 and test > 0:
#                 trials.append(
#                     [
#                         sampleTS,
#                         testTS,
#                         np.round(sampleTS / s1s, decimals=3),
#                         np.round(testTS / s1s, decimals=3),
#                         sample,
#                         test,
#                         rsps,
#                         np.round((testTS - sampleTS) / s1s) - 1,
#                     ]
#                 )
#                 sample = -1
#                 test = -1
#                 sampleTS = -1
#                 testTS = -1
#                 rsps = -1
#             lastCue = cue
#             lastTS = cueTS
# 
#         if (
#             test > 0
#             and rsps < 0
#             and events[eidx][0] >= testTS + s1s
#             and events[eidx][0] < lastTS + 2 * s1s
#             and (events[eidx][1] & 0x01) > 0
#         ):  # response window
#             rsps = 1
# 
#     if sample > 0 and test > 0:  # pick up last trial
#         trials.append(
#             [
#                 sampleTS,
#                 testTS,
#                 np.round(sampleTS / s1s, decimals=3),
#                 np.round(testTS / s1s, decimals=3),
#                 sample,
#                 test,
#                 rsps,
#                 np.round((testTS - sampleTS) / s1s) - 1,
#             ]
#         )
# 
#     return trials


def getEvents(path=""):
    syncs = np.load(os.path.join(path,"sync_raw.npy"))
    blockCount = 0
    ts = 0
    events = deque() 
    pct = 0
    while ts < (len(syncs)):
        currPct = ts * 100 // len(syncs)  # pct for tracking progress
        if currPct > pct + 9:
            print(currPct)
            pct = currPct

        if syncs[ts] == 0:
            blockCount += 1
            ts += 1
        else:
            if blockCount < 34:  # skip processed data
                ts += 1
            else:
                blockCount = 0  # triggerd new block
                state = [
                    ts,
                    0,
                    0,
                    0,
                    0,
                ]  # bit-masked behavior data from MCU, ref. 12F1572SyncEncoder
                state[1] = 1 if np.sum(syncs[ts + 5 : ts + 10]) > 128 else 0
                state[2] = 1 if np.sum(syncs[ts + 10 : ts + 14]) > 127 else 0
                state[3] = 1 if np.sum(syncs[ts + 14 : ts + 19]) > 128 else 0
                state[4] = 1 if np.sum(syncs[ts + 19 : ts + 24]) > 128 else 0
                events.append(state)
                ts += 24
    print(100)
    events = np.array(events)
    return events


def writeEvents(events, trials,path=""):
    np.save(os.path.join(path,"sync_events.npy"), events)
    np.save(os.path.join(path,"sync_trials.npy"), trials)


# def filter_events(events):
#     """
#     if consecutive identical events interrupted by minimal noise, clean up the noise
# 
#     Parameters
#     ----------
#     events : [TS, bit-masked-type]
#         as exported from neuropixels binary data
# 
#     Returns
#     -------
#     output : [TS, bit-masked-type]
#         cleaned events
#     """
#     if isinstance(events,list):
#         events=np.array(events)
# 
# 
#     lick_interval = 30000 * 0.05  # 50ms
#     cue_interval = 30000 * 0.5
#     prev_idx = np.argmax(events[:, 1])  # first non-zero
#     prev_TS = events[prev_idx, 0]
#     for i in range(prev_idx + 1, events.shape[0]):
#         if events[i, 1] == 1:
#             curr_TS = events[i, 0]
#             if curr_TS - prev_TS < lick_interval:
#                 events[prev_idx:i, 1] = 1
#             prev_idx = i
#             prev_TS = curr_TS
# 
#     prev_idx = np.argmax(events[:, 3])
#     prev_TS = events[prev_idx, 0]
#     for i in range(prev_idx + 1, events.shape[0]):
#         if events[i, 3] == 1:
#             curr_TS = events[i, 0]
#             if curr_TS - prev_TS < cue_interval:
#                 events[prev_idx:i, 3] = 1
#             prev_idx = i
#             prev_TS = curr_TS
# 
#     prev_idx = np.argmax(events[:, 4])
#     prev_TS = events[prev_idx, 0]
#     for i in range(prev_idx + 1, events.shape[0]):
#         if events[i, 4] == 1:
#             curr_TS = events[i, 0]
#             if curr_TS - prev_TS < cue_interval:
#                 events[prev_idx:i, 4] = 1
#             prev_idx = i
#             prev_TS = curr_TS
# 
#     events[:, 2] = 0
#     output = []
#     output.append([events[0, 0], events[0, 1] + events[0, 3] * 4 + events[0, 4] * 8])
#     for i in range(1, events.shape[0]):
#         val = events[i, 1] + events[i, 3] * 4 + events[i, 4] * 8
# 
#         if not val == output[-1][1]:
#             output.append([events[i, 0], val])
#     return output


def filter_imec(events,sig_type='cue',update_interval=45): # events:2-columns, TS and cue
    if sig_type not in ('cue','lick'):
        raise ValueError(f'Invalid signal type {sig_type}')
    buf=deque()
    buf.append(events[0,:])
    for ii in range(1,events.shape[0]):
        if events[ii,0]-events[ii-1,0]>update_interval:
            buf.append([events[ii-1,0]+30,0]);
        buf.append(events[ii,:])

    events=np.array(buf)

    edge_diff=np.diff(np.hstack(([0], events[:-1, 1], [0])))
    edge_on = events[np.where(edge_diff>0)[0],0]
    edge_off = events[np.where(edge_diff<0)[0],0]
    edge_sel = np.ones((edge_on.shape[0]), dtype=bool)
    if sig_type=='lick':
        edge_sel[(edge_off-edge_on)<0.01*30000]=False
        return (edge_on[edge_sel],edge_off[edge_sel])

    else:
        for idx in range(1, edge_on.shape[0] - 1):
            if (
                (edge_off[idx] - edge_on[idx]) / (edge_off[idx] - edge_off[idx - 1]) < 0.6
                and (edge_off[idx] - edge_on[idx]) / (edge_on[idx + 1] - edge_on[idx]) < 0.6
                and (edge_off[idx] - edge_on[idx]) < 0.02*30000
            ):
                edge_sel[idx] = False

        # first and last

        if ((edge_off[0] - edge_on[0]) / (edge_on[1] - edge_on[0]) < 0.6
            and (edge_off[0] - edge_on[0]) < 0.02*30000):
            edge_sel[0] = False

        if ((edge_off[-1] - edge_on[-1]) / (edge_on[-1] - edge_on[-2]) < 0.6
            and (edge_off[-1] - edge_on[-1]) < 0.02*30000):
            edge_sel[-1] = False
        edge_on = edge_on[edge_sel]
        edge_off = edge_off[edge_sel]

        return (edge_on,edge_off)

def events2trials(events,sync_type='zx'): # assume have been pre-processed
    if sync_type not in ('zx','xd'):
        raise ValueError(f'sync signal type error {sync_type}');
    if sync_type=='zx':
        lickon,lickoff=filter_imec(events[:,(0,1)],sig_type='lick',update_interval=90)
        e3on,e3off=filter_imec(events[:,(0,3)],sig_type='cue',update_interval=90)
        e4on,e4off=filter_imec(events[:,(0,4)],sig_type='cue',update_interval=90)
    else:
        lickon,lickoff=filter_imec(events[:,(0,1)],sig_type='lick',update_interval=45)
        e3on,e3off=filter_imec(events[:,(0,3)],sig_type='cue',update_interval=45)
        e4on,e4off=filter_imec(events[:,(0,4)],sig_type='cue',update_interval=45)

    e3on=np.vstack((e3on,np.full_like(e3on,fill_value=4)))
    e3off=np.vstack((e3off,np.full_like(e3off,fill_value=-4)))
    e4on=np.vstack((e4on,np.full_like(e4on,fill_value=8)))
    e4off=np.vstack((e4off,np.full_like(e4off,fill_value=-8)))
    evts=np.transpose(np.hstack((e3on,e3off,e4on,e4off)))

    sidx=np.lexsort((evts[:,1],evts[:,0]))
    evts=evts[sidx,:]

    s1s = 30000
    trials = deque()
    lastTS = -300000 * 20
    lastCue = -1
    cueTS = -1
    rsps = -1
    cue = -1
    sample = -1
    test = -1
    sampleTS = -1
    testTS = -1

    prev_4off=-1
    prev_8off=-1

    for eidx in range(evts.shape[0]):
        if evts[eidx,1]==-4:
            prev_4off=evts[eidx,1]
            continue
        elif evts[eidx,1]==-8:
            prev_8off=evts[eidx,1]
            continue
        elif evts[eidx,1]==4 and (evts[eidx,0]-prev_4off)<0.1*s1s:
            continue
        elif evts[eidx,1]==8 and (evts[eidx,0]-prev_8off)<0.1*s1s:
            continue
        cue=evts[eidx,1]
        cueTS = evts[eidx][0]  # Time stamp
        if cue > 0 and cueTS > lastTS + 1.1 * s1s and cueTS < lastTS + 2 * s1s:
            print("error processing evt idx ", eidx)
        elif (
            cue > 0 and cueTS > lastTS + s1s * 2 and cueTS < lastTS + s1s * 8
        ):  # following delay
            sample = lastCue
            sampleTS = lastTS
            test = cue
            testTS = cueTS

            lastCue = cue
            lastTS = cueTS

        elif cue > 0 and cueTS > lastTS + s1s * 8:  # following ITI
            if sample > 0 and test > 0:
                trials.append(
                    [
                        sampleTS,
                        testTS,
                        np.round(sampleTS / s1s, decimals=3),
                        np.round(testTS / s1s, decimals=3),
                        sample,
                        test,
                        rsps,
                        np.round((testTS - sampleTS) / s1s) - 1,
                    ]
                )
                sample = -1
                test = -1
                sampleTS = -1
                testTS = -1
                rsps = -1
            lastCue = cue
            lastTS = cueTS

    if sample > 0 and test > 0:  # pick up last trial
        trials.append(
            [
                sampleTS,
                testTS,
                np.round(sampleTS / s1s, decimals=3),
                np.round(testTS / s1s, decimals=3),
                sample,
                test,
                rsps,
                np.round((testTS - sampleTS) / s1s) - 1,
            ]
        )
    trials=np.array(trials)
    for ii in range(trials.shape[0]):
        if np.any((
                np.logical_and(lickon>trials[ii,1]+s1s,lickon<=trials[ii,1]+2*s1s),
                np.logical_and(lickoff>trials[ii,1]+s1s,lickoff<=trials[ii,1]+2*s1s),
                np.logical_and(lickon<trials[ii,1]+s1s,lickoff>=trials[ii,1]+2*s1s),
                )):
            trials[ii,6]=1

    return trials

def runsync():
    events = getEvents()
    events = filter_events(events)
    trials = parseDPAEvents(events)
    writeEvents(events, trials)
    return trials


def old_main():
    print(os.getcwd())
    breakpoint()
    if os.path.exists("sync_trials.npy"):
        print("data exist")
        sys.exit(0)

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        fl = glob.glob("*.ap.bin")
        if len(fl) != 1:
            print("Unexpected .ap.bin file condition")
            sys.exit(100)
        filename = fl[0]
    print(f"Processing file: {filename}")
    readsync(filename)

    events = getEvents()
    events = filter_events(events)
    trials = parseDPAEvents(events)
    writeEvents(events, trials)

# with open('events.csv','w',newline='\r\n') as

if __name__ == "__main__":  # Debuging entry. Import module when possible.
    if "snakemake" in globals():
        outf=snakemake.output[0]
        if outf.endswith('sync_trials.npy'):
            print(f"Processing file: {outf}")
            if 'LN2' in snakemake.wildcards[0]:
                events = np.load(os.path.join(snakemake.wildcards[0],'sync_events_5.npy')).astype(np.int32)
                trials = events2trials(events,sync_type='xd')
                writeEvents(events, trials,snakemake.wildcards[0])
            else:
                events = getEvents(snakemake.wildcards[0]).astype(np.int32)
                trials = events2trials(events,sync_type='zx')
                writeEvents(events, trials,snakemake.wildcards[0])
