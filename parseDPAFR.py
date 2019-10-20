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
        
def baselineVector(oneTS,tsId,trials):
    tIdices=range(trials.shape[0])
    base=[]
    for tIdx in tIdices:
        for binCount in np.histogram(oneTS[tsId==tIdx],bins=4,range=(-60000,-30000))[0]:
            base.append(binCount)
    
    if len(base)>0 and np.std(base):
        return (np.mean(base), np.std(base))
    else:
        print('Error calculating base vector unit#%d' % (tsId,))
        return (0,32767)

def toHist(trials,oneTS,tsId,sample,delay,baseVector=(0,1)):
    sel=np.nonzero(np.bitwise_and(trials[:,4]==sample , trials[:,7]==delay))[0]
    return (np.histogram(oneTS[np.isin(tsId,sel)],np.linspace(-60000,30000*(delay+6),num=(delay+8)*4+1))[0])/len(sel)

           


def alignHeatmap(spkTS,spkCluster,unitInfo,trials):
    heat43Raw=[]
    heat46Raw=[]
    heat83Raw=[]
    heat86Raw=[]
    
    baseVecAll=[]
    depth=[]
    
    spkNThresh=spkTS[-1]/s1s*2    
    
    for infoIdx in range(unitInfo.shape[0]):
        spkIdx=unitInfo['id'][infoIdx]
        wf=unitInfo.iloc[infoIdx,8]=='good' or (math.isnan(unitInfo.iloc[infoIdx,8]) and unitInfo.iloc[infoIdx,3]=='good')
        spkCount=unitInfo.iloc[infoIdx,9]
        if spkCount>spkNThresh and wf:
            oneTSAll=(spkTS[spkCluster==spkIdx]).astype('int64')
            (oneTS,tsId)=trialAlign(trials,oneTSAll)
            baseVec=baselineVector(oneTS,tsId,trials)
            baseVecAll.append(baseVec)

            heat43Raw.append(toHist(trials,oneTS,tsId,4,3))
            heat46Raw.append(toHist(trials,oneTS,tsId,4,6))
            heat83Raw.append(toHist(trials,oneTS,tsId,8,3))
            heat86Raw.append(toHist(trials,oneTS,tsId,8,6))

            depth.append(unitInfo['depth'][infoIdx])
            
    depth=np.array(depth)
    baseVecAll=np.array(baseVecAll)
    dIdx=np.argsort(depth)
     
    heat43Raw=np.array(heat43Raw)
    heat46Raw=np.array(heat46Raw)
    heat83Raw=np.array(heat83Raw)
    heat86Raw=np.array(heat86Raw)

    return ((heat43Raw[dIdx,:],heat46Raw[dIdx,:],heat83Raw[dIdx,:],heat86Raw[dIdx,:]),baseVecAll[dIdx],depth[dIdx])

def plotOne(data,delay,ax,ylbl):
    
    plt.imshow(data,cmap='jet',aspect='auto',vmin=-3,vmax=3)
    
    
    if delay==6:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in np.array([2,3,9,10])*4-0.5]
        ax.set_xticks(np.array([2,7,12])*4-0.5)
        ax.set_xticklabels([0,5,10])
#        ax.set_xlabel('Time (s)')
        
    elif delay==3:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in np.array([2,3,6,7])*4-0.5]
        ax.set_xticks(np.array([2,7])*4-0.5)
        ax.set_xticklabels([0,5])
    
    if ylbl:
        ax.set_ylabel('Unit #')


def plotOneSel(A,B,delay,ax,ylbl):
    
    plt.imshow((B-A)/(B+A),cmap='jet',aspect='auto',vmin=-1,vmax=1)

    if delay==6:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in np.array([2,3,9,10])*4-0.5]
        ax.set_xticks(np.array([2,7,10])*4-0.5)
        ax.set_xticklabels([0,5,10])
        
        
    elif delay==3:
        [plt.plot([x,x],ax.get_ylim(),'-w') for x in np.array([2,3,6,7])*4-0.5]
        ax.set_xticks(np.array([2,7])*4-0.5)
        ax.set_xticklabels([0,5])
    
    if ylbl:
        ax.set_ylabel('Unit #')    

    ax.set_xlabel('Time (s)')

def plotHeatmap(raw,base,depth):
    import os
    import re
    cwd=os.getcwd();
    grps=re.search('19.*(?=_cleaned)',cwd)

    
    fh=plt.figure(3,figsize=[7.5,10])
    ax=plt.subplot(3,3,1)
    plotOne(((raw[0].transpose()-base[:,0])/base[:,1]).transpose(),3,ax,True)
    ax.set_title('S1 3s delay')
    ax=plt.subplot(3,3,2)
    plotOne(((raw[2].transpose()-base[:,0])/base[:,1]).transpose(),3,ax,False)
    ax.set_title('S2 3s delay')
    ax=plt.subplot(3,3,4)
    plotOne(((raw[1].transpose()-base[:,0])/base[:,1]).transpose(),6,ax,True)
    ax.set_title('S1 6s delay')
    ax=plt.subplot(3,3,5)
    plotOne(((raw[3].transpose()-base[:,0])/base[:,1]).transpose(),6,ax,False)
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
#    import os
#    os.chdir('D:\Data\\191018-DPA-Learning5_28_g1\\191018-DPA-Learning5_28_g1_imec1_cleaned')
    
    s1s=30000
    spkTS=np.load("spike_times.npy")
    spkCluster=np.load("spike_clusters.npy")
    
    unitInfo=pd.read_csv('cluster_info.tsv',sep='\t')
    
    trials=np.empty([0])
    with h5py.File('events.hdf5','r') as fe:
        dset=fe['trials']
        trials=np.array(dset,dtype='int64')    
    (raw,baseVec,depth)=alignHeatmap(spkTS,spkCluster,unitInfo,trials)
    (fh,ax)=plotHeatmap(raw,baseVec,depth)