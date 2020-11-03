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
import scikits.bootstrap as boot
from mlxtend.evaluate import permutation_test
from matplotlib import rcParams
import pickle



def others():
    
    candidate_per_sess=[]
    
    
    
    key_count_congru={}
    key_count_incong={}
    reg_count={}
    reg_pattern_count={}
    # reg_count_congru={}
    #
    for sess_id in range(114):
        if not os.path.isfile(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id)):
            continue
        
        mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
      
        ## count candidiates
        
        pref=np.array(mat['prefs'])[:,0]
        rego=mat['regs']
        regs=[x[0][0] for x in rego]
        
        reg_set=np.unique(regs)
        reg_set=reg_set[reg_set!='Unlabeled']
     
        for one_reg in reg_set:
            if one_reg in reg_count:
                reg_count[one_reg]+=regs.count(one_reg)
            else:
                reg_count[one_reg]=regs.count(one_reg)
    
        
        filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
        if not os.path.isfile(filt_fpath):
            continue
    
    
    
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
    
    
        ## end of counting candidates
        
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
            #select only inter-region connection
            if np.unique(reg_grp).shape[0]<2:
                # print('local')
                continue
            if np.isin('Unlabeled',reg_grp):
                # print('Unlabeled')
                continue
            

    
            pref_grp=[pref[x] for x in patt['neurons']]
            if np.unique(pref_grp).shape[0]<2:
                # prefered_samp=np.unique(pref_grp)[0]
                if lagkey in key_count_congru:
                    key_count_congru[lagkey]+=1
                else:
                    key_count_congru[lagkey]=1
                    
                for one_reg in reg_grp:
                    if one_reg in reg_pattern_count:
                        reg_pattern_count[one_reg]+=1
                    else:
                        reg_pattern_count[one_reg]=1     
                
            else:
                # prefered_samp=0
                if lagkey in key_count_incong:
                    key_count_incong[lagkey]+=1
                else:
                    key_count_incong[lagkey]=1
    
    
    # exported from matlab
    
    (cong_candi,incong_candi)=(*np.sum(np.array(candidate_per_sess),axis=0)[1:3],)
    
    
    sys.exit()
    
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    
    (fh, ax) = plt.subplots(1, 1, figsize=(19 / 2.54, 4 / 2.54), dpi=300)
    
    for idx,key in enumerate(keys):
        mm=[key_count_incong[key]/incong_candi,key_count_congru[key]/cong_candi]
        ax.bar(idx-0.15+1,mm[0],width=0.3,color='k',edgecolor='k')
        ax.bar(idx+0.15+1,mm[1],width=0.3,color='w',edgecolor='k')
        # ax.errorbar([1,2],mm,np.vstack((incongru_boot,congru_boot)),color='none',ecolor='grey',capsize=3)
        # ax.set_yscale('log')
        # ax.set_ylim([1e-4,5e-2])
        # ax.set_ylabel('Pattern density')
        ax.set_xticks(np.arange(1,21))
        # ax.set_xticklabels(['incongru.','congruent'],rotation=45,ha='right')
    # plt.close('all')
    fh.savefig('4su_temproal_dist.pdf',bbox_inches='tight')
    
    
    
    (fh, ax) = plt.subplots(1, 1, figsize=(22 / 2.54, 4 / 2.54), dpi=300)
    idx=0
    reglabel=[]
    keyset=list(reg_count.keys())
    keyset.sort()
    for key in keyset:
        if reg_count[key]>=30:
            if key not in reg_pattern_count:
                reg_pattern_count[key]=0
            reglabel.append(key)
            ax.bar(idx,reg_pattern_count[key]/reg_count[key],width=0.7,color='w',edgecolor='k')
            idx+=1
    
    # ax.errorbar([1,2],mm,np.vstack((incongru_boot,congru_boot)),color='none',ecolor='grey',capsize=3)
    ax.set_yscale('log')
        # ax.set_ylim([1e-4,5e-2])
        # ax.set_ylabel('Pattern density')
    ax.set_xlim([0,idx])    
    ax.set_xticks(np.arange(len(reglabel))+1)
    ax.set_xticklabels(reglabel,rotation=90,ha='center',va='top')
    ax.set_ylabel('Patterns / neuron')
    # plt.close('all')
    fh.savefig('4su_spatial_dist.pdf',bbox_inches='tight')
    
    
            
    ####################
    
def tempo_hz(congru_stats):
    congru_per_key={}
    keys=[0,1,2,3,11,12,13,22,23,33,111,112,113,122,123,133,222,223,233,333]
    for key in keys:
        congru_per_key[key]={'data':[],'mean':None,'ci':None,'n':None,'sem':None}
    
    patt_keys=[np.int64(np.dot(x.magnitude//4,[100,10,1])) for x in congru_stats['lag']]
    
    for idx,patt_key in enumerate(patt_keys):
        congru_per_key[patt_key]['data'].append(np.mean(congru_stats['pertrial'][idx][0]/6))
    
    for key in keys:
        congru_per_key[key]['mean']=np.mean(congru_per_key[key]['data'])
        congru_per_key[key]['ci']=boot.ci(congru_per_key[key]['data'], np.mean,n_samples=1000)
        congru_per_key[key]['n']=len(congru_per_key[key]['data'])
        congru_per_key[key]['sem']=np.std(congru_per_key[key]['data'],ddof=1)/np.sqrt(len(congru_per_key[key]['data']))
        
    (fh, ax) = plt.subplots(1, 1, figsize=(10 / 2.54, 4 / 2.54), dpi=300)
    ax.bar(np.arange(20)+1,[x['mean'] for x in congru_per_key.values()],
           color='w',edgecolor='k')
    ax.errorbar(np.arange(20)+1,
                [x['mean'] for x in congru_per_key.values()],
                # np.array([x['ci'] for x in congru_per_key.values()]).transpose(),
                [x['sem'] for x in congru_per_key.values()],
                color='none',ecolor='k',capsize=3)
    ax.set_ylabel('Pattern freq. (Hz)')
    ax.set_xlabel('Pattern #')
    a2 = ax.twinx()
    a2.plot(np.arange(20)+1,[x['n'] for x in congru_per_key.values()],'-r')
    a2.set_ylabel('Pattern number',color='r')
    a2.set_ylim((0,1000))
    a2.set_yticks([0,500,1000])
    fh.savefig('4su_temproal_dist.pdf',bbox_inches='tight')
    
    
if __name__=='__main__':
    if 'congru_stats' not in dir():
        r=pickle.load(open('spade_stats.p','rb'))
        congru_stats=r['congru_stats']
    tempo_hz(congru_stats)
    