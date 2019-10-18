# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


#import glob
import h5py
import numpy as np
#import os




def parseGNGEvents(events):
    s1s=30000
    trials=[]
    lastTS=-300000
    lastCue=-1
    cueTS=-1
    rsps=-1
    cue=-1
    
    for eidx in range(len(events)):
        cue=events[eidx][1] & 0x0c
        cueTS=events[eidx][0]
        if cue>0 and cueTS>(lastTS+s1s):
            if lastCue>0:
                trials.append([lastTS,lastCue,rsps])
                rsps=-1;
            lastCue=cue;
            lastTS=cueTS;

        if (lastCue>0 and rsps<0 
            and events[eidx][0]>= lastTS+30000 
            and events[eidx][0]< lastTS+90000 
            and (events[eidx][1] & 0x03)>0):
            rsps=events[eidx][1] & 0x03
    
    return trials
        

def parseDPAEvents(events):
    """
    if no cue,  assume sample
    if has cue history check time diff
    time diff < 3 discard
    time diff > 3 < 10, this cue is test, last cue is sample
    time diff >10
    
    if sample and test set
        if responsed
            add responsed trial
        else
            add unresponsed trial
        
        reinit sample test response
        this cue is sample
"""
    s1s=30000
    trials=[]
    lastTS=-300000*20
    lastCue=-1
    cueTS=-1
    rsps=-1
    cue=-1
    sample=-1
    test=-1
    sampleTS=-1
    testTS=-1
    
    for eidx in range(len(events)):
        cue=events[eidx][1] & 0x0c
        cueTS=events[eidx][0]
        if cue>0 and cueTS>lastTS+s1s and cueTS<lastTS+2*s1s:
            print("error processing evt idx ",eidx)
        elif cue>0 and cueTS>lastTS+s1s*2 and cueTS<lastTS+s1s*8:
            sample=lastCue
            sampleTS=lastTS
            test=cue
            testTS=cueTS
            
            lastCue=cue
            lastTS=cueTS
    
        elif cue>0 and cueTS>lastTS+s1s*8:
            if sample > 0 and test >0:
                trials.append([sampleTS,\
                               testTS,\
                               np.round(sampleTS/s1s,decimals=3),\
                               np.round(testTS/s1s,decimals=3),\
                               sample,\
                               test,\
                               rsps,\
                               np.round((testTS-sampleTS)/s1s)-1])
                sample=-1
                test=-1
                sampleTS=-1
                testTS=-1
                rsps=-1
            lastCue=cue
            lastTS=cueTS
            

    
    
        if (test>0 and rsps<0 
            and events[eidx][0]>= testTS+s1s
            and events[eidx][0]< lastTS+2*s1s
            and (events[eidx][1] & 0x03)>0):
            rsps=1
        
        if sample > 0 and test >0:
            trials.append([sampleTS,\
                           testTS,\
                           np.round(sampleTS/s1s,decimals=3),\
                           np.round(testTS/s1s,decimals=3),\
                           sample,\
                           test,\
                           rsps,\
                           np.round((testTS-sampleTS)/s1s)-1])
        
        
    return trials





def getEvents():
    syncs=np.empty([0])
    with h5py.File('sync.hdf5','r') as fs:
        dset=fs['sync']
        syncs=np.array(dset,dtype='int8')    
        syncs=syncs[0]
    
    
    blockCount=0
    ts=0
    events=[]
    pct=0
    while ts<(len(syncs)):
        currPct=ts*100//len(syncs)
        if currPct>pct:
            print(currPct)
            pct=currPct
        
        if syncs[ts]==0:
            blockCount+=1
            ts+=1
        else:
            if blockCount<34:
                ts+=1
            else:
                blockCount=0
                state=0;
                state+=(1 if np.sum(syncs[ts+5:ts+10])>128 else 0)
                state+=(2 if np.sum(syncs[ts+10:ts+14])>127 else 0)
                state+=(4 if np.sum(syncs[ts+14:ts+19])>128 else 0)
                state+=(8 if np.sum(syncs[ts+19:ts+24])>128 else 0)
                if (not events) or events[-1][1]!=state:
                    events.append([ts,state])
                ts+=24
    return events

def writeEvents(events,trials):
    with h5py.File('events.hdf5','w') as fw:
        evtDset=fw.create_dataset('events',data=np.array(events,dtype='i1'))
        tDset=fw.create_dataset('trials',data=np.array(trials,dtype='i4'))
            


if __name__=="__main__":
    events=getEvents()            
    trials=parseDPAEvents(events)
    writeEvents(events,trials)
    

    


# with open('events.csv','w',newline='\r\n') as 
            
            
            
        
        
