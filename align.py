# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:37:14 2019

@author: Libra
"""
import csv
import os
import phylib.utils._misc as phyutil
import re
import numpy as np
import h5py

        
def traverse(path):
    for basepath, dirs, files in os.walk(path):
        if 'cluster_info.tsv' in files:
            yield basepath
            
def getTrackRegion(regionL,mid,date,imecNo):
    depthL=[]
    for row in regionL:
        if row[0]==date and row[1]==mid and row[2]==imecNo:
            depthL.append([row[3],int(row[4]),int(row[5])])
    return depthL
            

def matchDepth(depth,depthL):
    for row in depthL:
        if row[1] <= depth < row[2]:
#            print('match')
            return row[0]
    return 'Unlabeled'
            


def judgePerformance(trials):
    if trials.shape[0]>=40:
        correctResp=np.bitwise_xor(trials[:,4]==trials[:,5],trials[:,6]==1)
        inWindow=np.zeros((trials.shape[0],),dtype='bool')
        i=40;
        while i<trials.shape[0]:
            if np.sum(correctResp[i-40:i])>=32:
                inWindow[i-40:i]=1
            i+=1
        if np.sum(inWindow)>=40:  #Well Trained
            return ('wellTrained',3,trials[inWindow,:])
                   
        else:
            inWindow=np.zeros((trials.shape[0],),dtype='bool')
            licks=(trials[:,6]==1)
            i=40;
            while i<trials.shape[0]:
                if np.sum(licks[i-40:i])>=16:  #Learning
                    inWindow[i-40:i]=1
                i+=1
            if np.sum(inWindow)>=40:
                return('learning',2,trials[inWindow,:])
            elif np.sum(licks)<=trials.shape[0]//10:    #Passive
                return ('passive',0,trials)
            else:
                return('transition',1,trials) #Unlabled

 

           
regionL=[]
with open('K:/neupix/meta/Sites.csv',newline='') as cf:
    creader=csv.reader(cf,dialect='excel')
    for row in creader:
        regionL.append(row)

FR_Th=1.0        

miceIdA=re.compile('Learning\\d[-_](\d*)[_-]')
dateA=re.compile('(19\\d\\d\\d\\d)\\D')
imecNo=re.compile('imec(\\d)')

regionMatched=[]

for path in traverse('K:/neupix/DataSum/'):
    midGrps=miceIdA.search(path)
    dateGrps=dateA.search(path)
    imecGrps=imecNo.search(path)
    
    if midGrps and dateGrps and imecGrps:
        trials=np.empty([0])
        with h5py.File(os.path.join(path,'events.hdf5'),'r') as fe:
            dset=fe['trials']
            trials=np.array(dset,dtype='int32')
        (perfType,perfIdx,selTrials)=judgePerformance(trials)    
        
        depthL=getTrackRegion(regionL,midGrps.group(1),dateGrps.group(1),imecGrps.group(1))
        unitInfo=phyutil.read_tsv(os.path.join(path,'cluster_info.tsv'))
        spkTS=np.load(os.path.join(path,'spike_times.npy'))
        s1s=30000
        spkNThresh=spkTS[-1]/s1s*FR_Th  

        for row in unitInfo:
            if row['KSLabel']=='good' and row['n_spikes']>=spkNThresh:
                suDepth=row['depth']
                reg=matchDepth(suDepth,depthL)
                regionMatched.append([path,row['id'],suDepth,reg,perfIdx])
    else:
        print('error parsing neuron track record')
        
allRegions=[]
allPerfType=[]
for u in regionMatched:
    allRegions.append(u[3])
    allPerfType.append(u[4])
regionSet=list(set(allRegions))
count=np.zeros((len(regionSet),4))

for idx in range(len(allRegions)):
    count[regionSet.index(allRegions[idx]),allPerfType[idx]]+=1

for r in range(len(regionSet)):
    print('%s, %d, %d, %d' % (regionSet[r],count[r,0],count[r,2],count[r,3]))    
    
    
    
    

    
#    
#    os.chdir(path)
#    (goodU,mU)=countNeurons()
#    goodSum.append(goodU)
#    mUSum.append(mU)
#    paths.append(path)
#    return (goodSum,mUSum,paths)            