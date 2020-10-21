# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 10:34:34 2020

@author: Libra
"""
import numpy as np
import sys
import scipy.stats as stats
import scipy.io
import os.path
from scipy.special import comb
import matplotlib.pyplot as plt
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
    if not os.path.isfile(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id)):
        continue
    mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    pref=np.array(mat['prefs'])[:,0]
    rego=mat['regs']
    regs=[x[0][0] for x in rego]
    
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
    
    
    
    filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
    if not os.path.isfile(filt_fpath):
        continue

    trialInfo=mat['trialInfo']
    spktrials=mat['spiketrials']
    firingrate=mat['firingrate']
    
    s1sel=trialInfo[:,4]==4
    s2sel=trialInfo[:,4]==8
    
    goodsel=np.logical_and(trialInfo[:,8],trialInfo[:,9])
    errorsel=np.logical_not(trialInfo[:,9])
    
    
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

sys.exit()
        
###################

congru_dens=[]
incongru_dens=[]
for idx,cs in enumerate(candidate_per_sess):
    if cs[1]>1000:
        congru_dens.append(np.sum(np.array(congru_stats.sess_ids)==cs[0])/cs[1])
    if cs[2]>1000:
        incongru_dens.append(np.sum(np.array(incongru_stats.sess_ids)==cs[0])/cs[2])

import scikits.bootstrap as boot
from mlxtend.evaluate import permutation_test

congru_boot=boot.ci(congru_dens, np.mean,n_samples=1000)
incongru_boot=boot.ci(incongru_dens, np.mean,n_samples=1000)

p_value = permutation_test(congru_dens,incongru_dens,method='approximate',num_rounds=10000)


(fh, ax) = plt.subplots(1, 1, figsize=(2 / 2.54, 4 / 2.54), dpi=300)
mm=[np.mean(incongru_dens),np.mean(congru_dens)]
ax.bar(1,mm[0],color='k',edgecolor='k')
ax.bar(2,mm[1],color='w',edgecolor='k')
ax.errorbar([1,2],mm,np.vstack((incongru_boot,congru_boot)),color='none',ecolor='grey',capsize=3)
ax.set_yscale('log')
ax.set_ylim([1e-4,5e-2])
ax.set_ylabel('Pattern density')
ax.set_xticks([1,2])
ax.set_xticklabels(['incongru.','congruent'],rotation=45,ha='right')
# plt.close('all')
fh.savefig('spade_4su_pattern_density.pdf',bbox_inches='tight')


####################

s1sigsel=np.array(congru_stats.perHz_pvalues)[:,1]<0.05
s2sigsel=np.array(congru_stats.perHz_pvalues)[:,2]<0.05
(fh, ax) = plt.subplots(1, 1, figsize=(5 / 2.54, 5 / 2.54), dpi=300)
# error on yaxis
ax.scatter(np.array(congru_stats.perHz_mm)[np.logical_not(s1sigsel),1],
           np.array(congru_stats.perHz_mm)[np.logical_not(s1sigsel),0],
           s=1,c='silver',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.perHz_mm)[np.logical_not(s2sigsel),3],
           np.array(congru_stats.perHz_mm)[np.logical_not(s2sigsel),2],
           s=1,c='silver',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.perHz_mm)[s1sigsel,1],
           np.array(congru_stats.perHz_mm)[s1sigsel,0],
           s=1,c='r',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.perHz_mm)[s2sigsel,3],
           np.array(congru_stats.perHz_mm)[s2sigsel,2],
           s=1,c='r',marker='.',alpha=0.4)
ax.plot([0,0.26],[0,0.26],'--k')

ax.set_yticks([0,0.1,0.2])
ax.set_xticks([0,0.1,0.2])
ax.set_xlabel('patterns / spike / s, error trial')
ax.set_ylabel('patterns / spike / s, correct trial')
ax.set_xlim((0,0.26))
ax.set_ylim((0,0.26))
fh.savefig('spade_4su_pattern_correct_error.pdf',bbox_inches='tight')
#############################


s1sigsel=np.array(congru_stats.motif_pvalues)[:,1]<0.05
s2sigsel=np.array(congru_stats.motif_pvalues)[:,2]<0.05
(fh, ax) = plt.subplots(1, 1, figsize=(5 / 2.54, 5 / 2.54), dpi=300)
# error on yaxis
ax.scatter(np.array(congru_stats.mm)[np.logical_not(s1sigsel),1]/6,
           np.array(congru_stats.mm)[np.logical_not(s1sigsel),0]/6,
           s=1,c='silver',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.mm)[np.logical_not(s2sigsel),3]/6,
           np.array(congru_stats.mm)[np.logical_not(s2sigsel),2]/6,
           s=1,c='silver',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.mm)[s1sigsel,1]/6,
           np.array(congru_stats.mm)[s1sigsel,0]/6,
           s=1,c='r',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.mm)[s2sigsel,3]/6,
           np.array(congru_stats.mm)[s2sigsel,2]/6,
           s=1,c='r',marker='.',alpha=0.4)
ax.plot([0,3.6],[0,3.6],'--k')

ax.set_yticks(np.arange(0,4))
ax.set_xticks(np.arange(0,4))
ax.set_xlabel('patterns / s, error trial')
ax.set_ylabel('patterns / s, correct trial')
ax.set_xlim((0,3.6))
ax.set_ylim((0,3.6))
fh.savefig('spade_4su_raw_pattern_correct_error.pdf',bbox_inches='tight')












#############################
s1sel=np.array(congru_stats.prefered_samp)==1
s2sel=np.array(congru_stats.prefered_samp)==2
sigsel=np.array(congru_stats.perHz_pvalues)[:,0]<0.05
(fh, ax) = plt.subplots(1, 1, figsize=(5 / 2.54, 5 / 2.54), dpi=300)
#prefer 1, 1 on yaxis
ax.scatter(np.array(congru_stats.perHz_mm)[np.logical_and(s1sel, np.logical_not(sigsel)),2],
           np.array(congru_stats.perHz_mm)[np.logical_and(s1sel, np.logical_not(sigsel)),0],
           s=1,c='silver',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.perHz_mm)[np.logical_and(s2sel, np.logical_not(sigsel)),0],
           np.array(congru_stats.perHz_mm)[np.logical_and(s2sel, np.logical_not(sigsel)),2],
           s=1,c='silver',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.perHz_mm)[np.logical_and(s1sel, sigsel),2],
           np.array(congru_stats.perHz_mm)[np.logical_and(s1sel, sigsel),0],
           s=1,c='r',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.perHz_mm)[np.logical_and(s2sel, sigsel),0],
           np.array(congru_stats.perHz_mm)[np.logical_and(s2sel, sigsel),2],
           s=1,c='r',marker='.',alpha=0.4)
ax.plot([0,0.26],[0,0.26],'--k')

ax.set_yticks([0,0.1,0.2])
ax.set_xticks([0,0.1,0.2])
ax.set_xlabel('patterns / spike / s, non-prefered')
ax.set_ylabel('patterns / spike / s, prefered')
ax.set_xlim((0,0.26))
ax.set_ylim((0,0.26))
fh.savefig('spade_4su_pattern_prefered_nonprefered.pdf',bbox_inches='tight')



#####################################

s1sel=np.array(congru_stats.prefered_samp)==1
s2sel=np.array(congru_stats.prefered_samp)==2
sigsel=np.array(congru_stats.motif_pvalues)[:,0]<0.05
(fh, ax) = plt.subplots(1, 1, figsize=(5 / 2.54, 5 / 2.54), dpi=300)
#prefer 1, 1 on yaxis
ax.scatter(np.array(congru_stats.mm)[np.logical_and(s1sel, np.logical_not(sigsel)),2]/6,
           np.array(congru_stats.mm)[np.logical_and(s1sel, np.logical_not(sigsel)),0]/6,
           s=1,c='silver',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.mm)[np.logical_and(s2sel, np.logical_not(sigsel)),0]/6,
           np.array(congru_stats.mm)[np.logical_and(s2sel, np.logical_not(sigsel)),2]/6,
           s=1,c='silver',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.mm)[np.logical_and(s1sel, sigsel),2]/6,
           np.array(congru_stats.mm)[np.logical_and(s1sel, sigsel),0]/6,
           s=1,c='r',marker='.',alpha=0.4)

ax.scatter(np.array(congru_stats.mm)[np.logical_and(s2sel, sigsel),0]/6,
           np.array(congru_stats.mm)[np.logical_and(s2sel, sigsel),2]/6,
           s=1,c='r',marker='.',alpha=0.4)
ax.plot([0,3.6],[0,3.6],'--k')

ax.set_yticks(np.arange(0,4))
ax.set_xticks(np.arange(0,4))
ax.set_xlabel('patterns / s, non-prefered')
ax.set_ylabel('patterns / s, prefered')
ax.set_xlim((0,3.6))
ax.set_ylim((0,3.6))
fh.savefig('spade_4su_raw_pattern_prefered_nonprefered.pdf',bbox_inches='tight')


