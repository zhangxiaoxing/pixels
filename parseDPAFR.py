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

def trialAlign(trials, oneTS):
    oneTS=oneTS[np.bitwise_and(oneTS>=trials[0,0]-30000*5, oneTS<=(trials[-1,0]+30000*10))]
    TSidx=0
    tIdx=0
    tsId=np.ones_like(oneTS)
    while TSidx<len(oneTS) and tIdx<len(trials):
        if oneTS[TSidx]<trials[tIdx,0]+(trials[tIdx,7]+8)*30000:
            oneTS[TSidx]-=trials[tIdx,0]
            tsId[TSidx]=tIdx
            TSidx+=1
        else:
            tIdx+=1
    return (oneTS,tsId)
        
def baselineVector(oneTS,tsId):
    tIdices=np.unique(tsId)
    base=[]
    for tIdx in tIdices:
        base.append(np.sum(np.bitwise_and(np.bitwise_and(tsId==tIdx, oneTS>=-60000), oneTS<0))/4) #500ms bin
    if len(base)>0 and np.std(base):
        return (np.mean(base), np.std(base))
    else:
        return (0,32767)

def toHist(trials,oneTS,tsId,sample,delay,baseVector=(0,1)):
    sel=np.nonzero(np.bitwise_and(trials[:,4]==sample , trials[:,7]==delay))[0]
    return ((np.histogram(oneTS[np.isin(tsId,sel)],np.linspace(-60000,30000*(delay+6),num=(delay+8)*2+1))[0])/len(sel)-baseVector[0])/baseVector[1]

           


def alignHeatmap(spkTS,spkCluster,unitInfo,trials):
    heat43Raw=[]
    heat46Raw=[]
    heat83Raw=[]
    heat86Raw=[]
    
    heat43=[]
    heat46=[]
    heat83=[]
    heat86=[]
    depth=[]
    
    spkNThresh=spkTS[-1]/s1s*2    
    
    for infoIdx in range(unitInfo.shape[0]):
        spkIdx=unitInfo['id'][infoIdx]
        wf=unitInfo.iloc[infoIdx,8]=='good' or (math.isnan(unitInfo.iloc[infoIdx,8]) and unitInfo.iloc[infoIdx,3]=='good')
        spkCount=unitInfo.iloc[infoIdx,9]
        if spkCount>spkNThresh and wf:
            oneTSAll=(spkTS[spkCluster==spkIdx]).astype('int64')
            (oneTS,tsId)=trialAlign(trials,oneTSAll)
            baseVec=baselineVector(oneTS,tsId)
#            sel=np.nonzero(np.bitwise_and(trials[:,4]==4 , trials[:,7]==3))[0]
            heat43Raw.append(toHist(trials,oneTS,tsId,4,3))
            heat46Raw.append(toHist(trials,oneTS,tsId,4,6))
            heat83Raw.append(toHist(trials,oneTS,tsId,8,3))
            heat86Raw.append(toHist(trials,oneTS,tsId,8,6))


            heat43.append(toHist(trials,oneTS,tsId,4,3,baseVec))
            heat46.append(toHist(trials,oneTS,tsId,4,6,baseVec))
            heat83.append(toHist(trials,oneTS,tsId,8,3,baseVec))
            heat86.append(toHist(trials,oneTS,tsId,8,6,baseVec))
            depth.append(unitInfo['depth'][infoIdx])
            
    depth=np.array(depth)
    dIdx=np.argsort(depth)
    
    heat43=np.array(heat43)
    heat46=np.array(heat46)
    heat83=np.array(heat83)
    heat86=np.array(heat86)
            
    heat43Raw=np.array(heat43Raw)
    heat46Raw=np.array(heat46Raw)
    heat83Raw=np.array(heat83Raw)
    heat86Raw=np.array(heat86Raw)


    return ((heat43Raw[dIdx,:],heat46Raw[dIdx,:],heat83Raw[dIdx,:],heat86Raw[dIdx,:]),(heat43[dIdx,:],heat46[dIdx,:],heat83[dIdx,:],heat86[dIdx,:]),depth[dIdx])

def plotOne(data,delay,ax,ylbl):
    
    plt.imshow(data,cmap='jet',aspect='auto',vmin=-3,vmax=3)
    
    
    if delay==6:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in [3.5,5.5,17.5,19.5]]
        ax.set_xticks([3.5,13.5,23.5])
        ax.set_xticklabels([0,5,10])
#        ax.set_xlabel('Time (s)')
        
    elif delay==3:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in [3.5,5.5,11.5,13.5]]
        ax.set_xticks([3.5,13.5])
        ax.set_xticklabels([0,5])
    
    if ylbl:
        ax.set_ylabel('Unit #')


def plotOneSel(A,B,delay,ax,ylbl):
    
    plt.imshow((B-A)/(B+A),cmap='jet',aspect='auto',vmin=-1,vmax=1)

    if delay==6:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in [3.5,5.5,17.5,19.5]]
        ax.set_xticks([3.5,13.5,23.5])
        ax.set_xticklabels([0,5,10])
        
        
    elif delay==3:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in [3.5,5.5,11.5,13.5]]
        ax.set_xticks([3.5,13.5])
        ax.set_xticklabels([0,5])
    
    if ylbl:
        ax.set_ylabel('Unit #')    

    ax.set_xlabel('Time (s)')

def plotHeatmap(raw,normed,depth):
    import os
    import re
    cwd=os.getcwd();
    grps=re.search('19.*(?=_cleaned)',cwd)

    
    fh=plt.figure(3,figsize=[7.5,10])
    ax=plt.subplot(3,3,1)
    plotOne(normed[0],3,ax,True)
    ax.set_title('S1 3s delay')
    ax=plt.subplot(3,3,2)
    plotOne(normed[2],3,ax,False)
    ax.set_title('S2 3s delay')
    ax=plt.subplot(3,3,4)
    plotOne(normed[1],6,ax,True)
    ax.set_title('S1 6s delay')
    ax=plt.subplot(3,3,5)
    plotOne(normed[3],6,ax,False)
    ax.set_title('S2 6s delay')
    #depth plot
    ax=plt.subplot(1,3,3)
    plt.plot(3840-depth)
    ax.set_ylabel('depth (um)')
    ax.set_xlabel('unit #')
    plt.minorticks_on();
    plt.grid(b=True,which='both')
    
    #selectivity
    ax=plt.subplot(3,3,7)
    plotOneSel(raw[0],raw[2],3,ax,True)
    ax.set_title('3s selectivity')    

    ax=plt.subplot(3,3,8)
    plotOneSel(raw[1],raw[3],6,ax,False)
    ax.set_title('6s selectivity')        
    
    fh.suptitle(grps.group().replace('_cleaned',''))
    plt.tight_layout(rect=[0,0,1,0.95])
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
    (raw,normed,depth)=alignHeatmap(spkTS,spkCluster,unitInfo,trials)
    (fh,ax)=plotHeatmap(raw,normed,depth)