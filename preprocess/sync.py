# -*- coding: utf-8 -*-
"""
@author: Libra
"""

# import glob
import sys, os, glob
import numpy as np
from readsync import readsync


def parseGNGEvents(events):
    '''
    for Go/Nogo task event parsing
    '''
    s1s = 30000 #sample per seconds
    trials = []
    lastTS = -300000
    lastCue = -1
    cueTS = -1
    rsps = -1
    cue = -1

    for eidx in range(len(events)):
        cue = events[eidx][1] & 0x0C #bit mask S1, S2
        cueTS = events[eidx][0] #Time stamp
        if cue > 0 and cueTS > (lastTS + s1s): #update trial
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
                and (events[eidx][1] & 0x03) > 0 #bit mask lick L/R
        ):
            rsps = events[eidx][1] & 0x03

    return trials


def parseDPAEvents(events):
    """
    for DPA task parsing
    if no cue,  assume sample
    if has cue history check time diff
    time diff < 3 discard
    time diff > 3 < 10, assume current cue as test, previous cue as sample
    time diff >10

    if sample and test assigned
        if responsed
            add responsed trial
        else
            add unresponsed trial

        reinit sample test response
        assume current cue as sample
"""
    s1s = 30000
    trials = []
    lastTS = -300000 * 20
    lastCue = -1
    cueTS = -1
    rsps = -1
    cue = -1
    sample = -1
    test = -1
    sampleTS = -1
    testTS = -1

    for eidx in range(len(events)):
        cue = events[eidx][1] & 0x0C #bit mask S1, S2
        cueTS = events[eidx][0] #Time stamp
        if cue > 0 and cueTS > lastTS + 1.1 * s1s and cueTS < lastTS + 2 * s1s:
            print("error processing evt idx ", eidx)
        elif cue > 0 and cueTS > lastTS + s1s * 2 and cueTS < lastTS + s1s * 8: # following delay
            sample = lastCue
            sampleTS = lastTS
            test = cue
            testTS = cueTS

            lastCue = cue
            lastTS = cueTS

        elif cue > 0 and cueTS > lastTS + s1s * 8: #following ITI
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

        if (test > 0 and rsps < 0 and events[eidx][0] >= testTS + s1s and events[eidx][0] < lastTS + 2 * s1s and (
                events[eidx][1] & 0x01) > 0): #response window
            rsps = 1

    if sample > 0 and test > 0: #pick up last trial
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

    return trials


def getEvents():
    syncs = np.load('sync_raw.npy') 
    blockCount = 0
    ts = 0
    events = []
    pct = 0
    while ts < (len(syncs)):
        currPct = ts * 100 // len(syncs) # pct for tracking progress
        if currPct > pct + 9:
            print(currPct)
            pct = currPct

        if syncs[ts] == 0:
            blockCount += 1
            ts += 1
        else:
            if blockCount < 34: #skip processed data
                ts += 1
            else:
                blockCount = 0 # triggerd new block
                state = [ts, 0, 0, 0, 0] # bit-masked behavior data from MCU, ref. 12F1572SyncEncoder
                state[1] = 1 if np.sum(syncs[ts + 5: ts + 10]) > 128 else 0
                state[2] = 1 if np.sum(syncs[ts + 10: ts + 14]) > 127 else 0
                state[3] = 1 if np.sum(syncs[ts + 14: ts + 19]) > 128 else 0
                state[4] = 1 if np.sum(syncs[ts + 19: ts + 24]) > 128 else 0
                if (not events) or not np.array_equal(events[-1][1:], state[1:]): #skip contineous event
                    events.append(state)
                ts += 24
    events = np.array(events)
    return events


def writeEvents(events, trials):
    np.save("sync_events.npy", events)
    np.save("sync_trials.npy", trials)
#    with h5py.File("events.hdf5", "w") as fw:
#        evtDset = fw.create_dataset("events", data=np.array(events, dtype="i4"))
#        tDset = fw.create_dataset("trials", data=np.array(trials, dtype="i4"))


def filter_events(events):
    '''
    if consecutive identical events interrupted by minimal noise, clean up the noise

    Parameters
    ----------
    events : [TS, bit-masked-type]
        as exported from neuropixels binary data

    Returns
    -------
    output : [TS, bit-masked-type]
        cleaned events

    '''
    lick_interval = 30000 * 0.05  # 50ms
    cue_interval = 30000 * 0.5
    prev_idx = np.argmax(events[:, 1])
    prev_TS = events[prev_idx, 0]
    for i in range(prev_idx + 1, events.shape[0]):
        if events[i, 1] == 1:
            curr_TS = events[i, 0]
            if curr_TS - prev_TS < lick_interval:
                events[prev_idx:i, 1] = 1
            prev_idx = i
            prev_TS = curr_TS

    prev_idx = np.argmax(events[:, 3])
    prev_TS = events[prev_idx, 0]
    for i in range(prev_idx + 1, events.shape[0]):
        if events[i, 3] == 1:
            curr_TS = events[i, 0]
            if curr_TS - prev_TS < cue_interval:
                events[prev_idx:i, 3] = 1
            prev_idx = i
            prev_TS = curr_TS

    prev_idx = np.argmax(events[:, 4])
    prev_TS = events[prev_idx, 0]
    for i in range(prev_idx + 1, events.shape[0]):
        if events[i, 4] == 1:
            curr_TS = events[i, 0]
            if curr_TS - prev_TS < cue_interval:
                events[prev_idx:i, 4] = 1
            prev_idx = i
            prev_TS = curr_TS
    events[:, 2] = 0
    output = []
    output.append([events[0, 0], events[0, 1] + events[0, 3] * 4 + events[0, 4] * 8])
    for i in range(1, events.shape[0]):
        val = events[i, 1] + events[i, 3] * 4 + events[i, 4] * 8
        if not val == output[-1][1]:
            output.append([events[i, 0], val])
    return output


def runsync():
    events = getEvents()
    events = filter_events(events)
    trials = parseDPAEvents(events)
    writeEvents(events, trials)
    return trials


if __name__ == "__main__": # Debuging entry. Import module when possible.
    if os.path.exists('sync_trials.npy'):
        print("data exist")
        sys.exit(0)

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        fl=glob.glob("*.ap.bin")
        if len(fl)!= 1:
            print("Unexpected .ap.bin file condition")
            sys.exit(100);
        filename = fl[0]
    print(f"Processing file: {filename}")
    readsync(filename)

    events = getEvents()
    events = filter_events(events)
    trials = parseDPAEvents(events)
    writeEvents(events, trials)

# with open('events.csv','w',newline='\r\n') as
