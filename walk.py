# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 17:12:15 2019

@author: Libra
"""

import os
#import parseDPAFR as pdpa
import phylib.utils._misc as phyutil
import numpy as np
import h5py
FR_Th=1.0

def traverse(path):
    for basepath, dirs, files in os.walk(path):
        if 'events.hdf5' in files:
            yield basepath
            
def countNeurons():
    unitInfo=phyutil.read_tsv('cluster_info.tsv')
    spkTS=np.load("spike_times.npy")
    s1s=30000
    spkNThresh=spkTS[-1]/s1s*FR_Th
    goodU=0
    mU=0
    for infoIdx in range(len(unitInfo)):
        wf=unitInfo[infoIdx].get('group')=='good' or ((unitInfo[infoIdx].get('group') is None) and unitInfo[infoIdx].get('KSLabel')=='good')
        spkCount=unitInfo[infoIdx]['n_spikes']
        if spkCount>spkNThresh and wf:
            goodU+=1
        elif spkCount>spkNThresh:
            mU+=1
    return (goodU,mU)

def judgePerformance(trials):
    if trials.shape[0]>=40:
        correctResp=np.logical_xor(trials[:,4]==trials[:,5],trials[:,6]==1)
        inWindow=np.zeros((trials.shape[0],),dtype='bool')
        i=40;
        while i<trials.shape[0]:
            if np.sum(correctResp[i-40:i])>=32:
                inWindow[i-40:i]=1
            i+=1
        if np.sum(inWindow)>=40:  #Well Trained
            return ('wellTrained',trials[inWindow,:])
                   
        else:
            inWindow=np.zeros((trials.shape[0],),dtype='bool')
            licks=(trials[:,6]==1)
            i=40;
            while i<trials.shape[0]:
                if np.sum(licks[i-40:i])>=16:  #Learning
                    inWindow[i-40:i]=1
                i+=1
            if np.sum(inWindow)>=40:
                return('learning',trials[inWindow,:])
            elif np.sum(licks)<=trials.shape[0]//10:    #Passive
                return ('passive',trials)
            else:
                return('transition',trials) #Unlabled
                
        
    
        
def walkNeurons():
    goodSum=[]
    mUSum=[]
    paths=[]
    perfTypes=[]
    for path in traverse('K:/neupix/DataSum/'):
        os.chdir(path)
        trials=np.empty([0])
        with h5py.File('events.hdf5','r') as fe:
            dset=fe['trials']
            trials=np.array(dset,dtype='int32')
        (perfType,selTrials)=judgePerformance(trials)
        (goodU,mU)=countNeurons()
        goodSum.append(goodU)
        mUSum.append(mU)
        paths.append(path)
        perfTypes.append(perfType)
    return (goodSum,mUSum,paths,perfTypes)

if __name__=="__main__":
    walkNeurons()    


#for path in traverse('K:/neupix/DataSum/'):
#    os.chdir(path)
    
