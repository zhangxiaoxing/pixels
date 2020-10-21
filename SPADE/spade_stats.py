# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 10:34:34 2020

@author: Libra
"""
import numpy as np
import scipy.stats as stats
import scipy.io
import os.path
from scipy.special import comb
# import neo

class results:
    def __init__(self):
        self.mm=[]
        self.perHz_mm=[]
        self.stdd=[]
        self.perHz_std=[]
        self.motif_pvalues=[]
        self.perHz_pvalues=[]
        # self.spikecount=[]
        self.sess_ids=[]
        self.prefered_samp=[]
        self.neu=[]
        self.neu_regs=[]


candidate_per_sess=[]


congru_stats=results();
incongru_stats=results();


for sess_id in range(114):
    
    filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
    if not os.path.isfile(filt_fpath):
        continue

    mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    trialInfo=mat['trialInfo']
    spktrials=mat['spiketrials']
    firingrate=mat['firingrate']
    pref=np.array(mat['prefs'])[:,0]
    rego=mat['regs']
    regs=[x[0][0] for x in rego]
    s1sel=trialInfo[:,4]==4
    s2sel=trialInfo[:,4]==8
    
    goodsel=np.logical_and(trialInfo[:,8],trialInfo[:,9])
    errorsel=np.logical_not(trialInfo[:,9])
    
    reg_set=np.unique(regs)
    
    reg_set=reg_set[reg_set!='Unlabeled']
    
    local_congru=0
    local_incongru=0
    
    total_n=comb(np.sum([x!='Unlabeled' for x in regs]),4)
    total_congru1=comb(np.sum(np.logical_and(pref==1,[x!='Unlabeled' for x in regs])),4)
    total_congru2=comb(np.sum(np.logical_and(pref==2,[x!='Unlabeled' for x in regs])),4)
    local_congru=0
    local_incongru=0
    for one_reg in reg_set:
        n_su=np.sum([x==one_reg for x in regs])
        local1=np.sum(np.logical_and(pref==1, [x==one_reg for x in regs]))
        local2=np.sum(np.logical_and(pref==2, [x==one_reg for x in regs]))
        local_congru+=(comb(local1,4)+comb(local2,4))
        local_incongru+=(comb(n_su,4)-comb(local1,4)-comb(local2,4))
        
        
    
    candidate_per_sess.append([sess_id,
                               total_congru1+total_congru2-local_congru,
                               total_n-total_congru1-local_incongru])
    
    
    
    
    
    
    r=np.load(filt_fpath,allow_pickle=True)
    patterns=r[0]
    
    for idx, patt in enumerate(patterns):
        reg_grp=[regs[x] for x in patt['neurons']]
        if np.unique(reg_grp).shape[0]<2:
            # print('local')
            continue
        if np.isin('Unlabeled',reg_grp):
            # print('Unlabeled')
            continue
        
        pref_grp=[pref[x] for x in patt['neurons']]
        if np.unique(pref_grp).shape[0]<2:
            currstat=congru_stats
            prefered_samp=np.unique(pref_grp)[0]
        else:
            currstat=incongru_stats
            prefered_samp=0
            
        currstat.sess_ids.append(sess_id)
        currstat.neu.append(patt['neurons'])
        currstat.prefered_samp.append(prefered_samp)
            
        currstat.neu_regs.append(reg_grp)
        times=patt['times']
        hist,edge=np.histogram(times,np.arange(0,7000*(len(trialInfo)+1),7000))
        # mm(idx)=(np.mean(hist[0:60]),np.mean(hist[60:120]))
        # tt=stats.ttest_ind(hist[0:60],hist[60:120])
        pv=[1,1,1]
        wrs=stats.ranksums(hist[np.logical_and(s1sel,goodsel)],hist[np.logical_and(s2sel,goodsel)])
        pv[0]=wrs[1]
        wrs=stats.ranksums(hist[np.logical_and(s1sel,goodsel)],hist[np.logical_and(s1sel,errorsel)])
        pv[1]=wrs[1]
        wrs=stats.ranksums(hist[np.logical_and(s2sel,goodsel)],hist[np.logical_and(s2sel,errorsel)])
        pv[2]=wrs[1]
        currstat.motif_pvalues.append(pv)
        
        
        
        currstat.mm.append([np.mean(hist[np.logical_and(s1sel,goodsel)]),
                   np.mean(hist[np.logical_and(s1sel,errorsel)]),
                   np.mean(hist[np.logical_and(s2sel,goodsel)]),
                   np.mean(hist[np.logical_and(s2sel,errorsel)])
                   ])
        
        currstat.stdd.append([np.std(hist[np.logical_and(s1sel,goodsel)]),
               np.std(hist[np.logical_and(s1sel,errorsel)]),
               np.std(hist[np.logical_and(s2sel,goodsel)]),
               np.std(hist[np.logical_and(s2sel,errorsel)])
               ])
        
        #firing rate
        
        
        s1g=np.divide(hist[np.logical_and(s1sel,goodsel)],[np.maximum(0.01,stats.gmean([firingrate[onesu][trl] for onesu in patt['neurons']])) for trl in np.nonzero(np.logical_and(s1sel,goodsel))[0]])
        s1e=np.divide(hist[np.logical_and(s1sel,errorsel)],[np.maximum(0.01,stats.gmean([firingrate[onesu][trl] for onesu in patt['neurons']])) for trl in np.nonzero(np.logical_and(s1sel,errorsel))[0]])
        s2g=np.divide(hist[np.logical_and(s2sel,goodsel)],[np.maximum(0.01,stats.gmean([firingrate[onesu][trl] for onesu in patt['neurons']])) for trl in np.nonzero(np.logical_and(s2sel,goodsel))[0]])
        s2e=np.divide(hist[np.logical_and(s2sel,errorsel)],[np.maximum(0.01,stats.gmean([firingrate[onesu][trl] for onesu in patt['neurons']])) for trl in np.nonzero(np.logical_and(s2sel,errorsel))[0]])
        
        currstat.perHz_mm.append([np.mean(s1g),
                         np.mean(s1e),
                         np.mean(s2g),
                         np.mean(s2e)])
        
        
        currstat.perHz_std.append([np.std(s1g),
                 np.std(s1e),
                 np.std(s2g),
                 np.std(s2e)])
        
            
        pv=[1,1,1]
        wrs=stats.ranksums(s1g,s2g)
        pv[0]=wrs[1]
        wrs=stats.ranksums(s1g,s1e)
        pv[1]=wrs[1]
        wrs=stats.ranksums(s2g,s2e)
        pv[2]=wrs[1]
        currstat.perHz_pvalues.append(pv)
        
        
        
        
        
        
# np.savez('spade_4su_stats.npz',congru_stats=congru_stats,incongru_stats=incongru_stats,allow_pickle=True)
        
        
# pv_arr=np.array(motif_pvalues)
        
        
    # duo_sel_set=np.nonzero(np.logical_and(motif_pvalues[:,0]<0.05,np.logical_or(motif_pvalues[:,1]<0.05,motif_pvalues[:,2]<0.05)))
    
    # selec_idces=np.nonzero(motif_pvalues<0.05)

# spkr=neo.io.NeoMatlabIO(filename='spkt_23_120_6.mat')   
# bl=spkr.read_block()                                            
# data=bl.segments[0].spiketrains                              




