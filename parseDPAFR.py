# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 17:24:12 2019

@author: Libra
"""

import numpy as np
import pandas as pd
import h5py
import math
import matplotlib.pyplot as plt

def trialAlign(trials,oneTS):
    oneTS=oneTS[np.bitwise_and(oneTS>=trials[0,0]-30000*5, oneTS<=(trials[-1,0]+30000*10))]
    TSidx=0
    tIdx=0;
    tsInfo=np.ones_like(oneTS);
    while TSidx<len(oneTS) and tIdx<len(trials):
        if oneTS[TSidx]<trials[tIdx,0]+(trials[tIdx,7]+8)*30000:
            oneTS[TSidx]-=trials[tIdx,0]
            tsInfo[TSidx]=tIdx
            TSidx+=1
        else:
            tIdx+=1
    return (oneTS,tsInfo)
        
def baselineVector(oneTS,tsInfo):
    tIdices=np.unique(tsInfo)
    base=[]
    for tIdx in tIdices:
        base.append(np.sum(np.bitwise_and(np.bitwise_and(tsInfo==tIdx, oneTS>=-60000), oneTS<0))/4) #500ms bin
    if len(base)>0 and np.std(base):
        return (np.mean(base), np.std(base))
    else:
        return (0,32767)

def toHist(trials,oneTS,tsInfo,sample,delay,baseVector):
    sel=np.nonzero(np.bitwise_and(trials[:,4]==sample , trials[:,7]==delay))[0]
    return ((np.histogram(oneTS[np.isin(tsInfo,sel)],np.linspace(-60000,30000*(delay+6),num=(delay+8)*2+1))[0])/len(sel)-baseVector[0])/baseVector[1]

           


def alignHeatmap(spkTS,spkCluster,unitInfo,trials):
    heat43=[]
    heat46=[]
    heat83=[]
    heat86=[]
    
    spkNThresh=spkTS[-1]/s1s*2    
    
    for idx in range(unitInfo.shape[0]):
        wf=unitInfo.iloc[idx,8]=='good' or (math.isnan(unitInfo.iloc[idx,8]) and unitInfo.iloc[idx,3]=='good')
        spkCount=unitInfo.iloc[idx,9]
        if spkCount>spkNThresh and wf:
            oneTSAll=(spkTS[spkCluster==idx]).astype('int64')
            (oneTS,tsInfo)=trialAlign(trials,oneTSAll)
            baseVec=baselineVector(oneTS,tsInfo)
#            sel=np.nonzero(np.bitwise_and(trials[:,4]==4 , trials[:,7]==3))[0]
            heat43.append(toHist(trials,oneTS,tsInfo,4,3,baseVec))
            heat46.append(toHist(trials,oneTS,tsInfo,4,6,baseVec))
            heat83.append(toHist(trials,oneTS,tsInfo,8,3,baseVec))
            heat86.append(toHist(trials,oneTS,tsInfo,8,6,baseVec))
    return (heat43,heat46,heat83,heat86)

def plotOne(data,delay,ax,ylbl):
    
    plt.imshow(data,cmap='jet',aspect='auto',vmin=-3,vmax=3)
    
    
    if delay==6:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in [3.5,5.5,17.5,19.5]]
        ax.set_xticks([3.5,13.5,23.5])
        ax.set_xticklabels([0,5,10])
        ax.set_xlabel('Time (s)')
        
    elif delay==3:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in [3.5,5.5,11.5,13.5]]
        ax.set_xticks([3.5,13.5])
        ax.set_xticklabels([0,5])
    
    if ylbl:
        ax.set_ylabel('Unit #')
    

def plotHeatmap(heat43,heat46,heat83,heat86):
#    plt.figure(0,figsize=[15,15])
#    plt.imshow(np.array(heat43),cmap='jet',vmin=-3,vmax=3)
#    plt.figure(1,figsize=[15,15])
#    plt.imshow(np.array(heat83),cmap='jet',vmin=-3,vmax=3)         
#    plt.figure(2,figsize=[15,15])
#    plt.imshow(np.array(heat46),cmap='jet',vmin=-3,vmax=3)
    fh=plt.figure(3,figsize=[7.5,10])
    ax=plt.subplot(2,2,1)
    plotOne(np.array(heat43),3,ax,True);
    ax=plt.subplot(2,2,2)
    plotOne(np.array(heat83),3,ax,False);
    ax=plt.subplot(2,2,3)
    plotOne(np.array(heat46),6,ax,True);
    ax=plt.subplot(2,2,4)
    plotOne(np.array(heat86),6,ax,False);
    plt.show();
    fh.savefig('heatmap.png',dpi=300,bbox_inches='tight')
    return (fh,ax)



if __name__=="__main__":
    s1s=30000
    spkTS=np.load("spike_times.npy")
    spkCluster=np.load("spike_clusters.npy")
    
    unitInfo=pd.read_csv('cluster_info.tsv',sep='\t')
    
    trials=np.empty([0])
    with h5py.File('events.hdf5','r') as fe:
        dset=fe['trials']
        trials=np.array(dset,dtype='int64')    
    (heat43,heat46,heat83,heat86)=alignHeatmap(spkTS,spkCluster,unitInfo,trials)
    (fh,ax)=plotHeatmap(heat43,heat46,heat83,heat86)