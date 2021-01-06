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



def perPattRatio():
    
    key_count_congru={}
    key_count_incong={}
    keys=[0,1,2,3,11,12,13,22,23,33,111,112,113,122,123,133,222,223,233,333]
    for key in keys:
        key_count_congru[key]=[]
        key_count_incong[key]=[]
    
    
    for sess_id in range(1,114):
        sess_congru={}
        sess_incong={}
        for key in keys:
            sess_congru[key]=0
            sess_incong[key]=0
        if not os.path.isfile(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id)):
            print(f'missing sessid {sess_id}')
            continue
    
        filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
        if not os.path.isfile(filt_fpath):
            continue
    
        mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    
        ## count candidiates
    
        pref=np.array(mat['prefs'])[:,0]
        rego=mat['regs']
        regs=[x[0][0] for x in rego]
    
        reg_set=np.unique(regs)
        reg_set=reg_set[reg_set!='Unlabeled']
    
    
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
    
            lagkey=(patt['lags'].magnitude[0]/4*100
                    +patt['lags'].magnitude[1]/4*10
                    +patt['lags'].magnitude[2]/4)
    
            pref_grp=[pref[x] for x in patt['neurons']]
            if np.unique(pref_grp).shape[0]<2:
                sess_congru[lagkey]+=1
    
            else:
                sess_incong[lagkey]+=1
        
        for key in keys:
            key_count_congru[key].append(sess_congru[key])
            key_count_incong[key].append(sess_incong[key])
    
    # exported from matlab
    
    
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    
    (fh, ax) = plt.subplots(1, 1, figsize=(6 / 2.54, 4 / 2.54), dpi=300)
    
    for idx,key in enumerate(keys):
        mm=[np.mean(key_count_incong[key]),np.mean(key_count_congru[key])]
        
        ax.bar(idx-0.2+1,mm[0],width=0.35,color='k',edgecolor='k',lw=0.5)
        ax.bar(idx+0.2+1,mm[1],width=0.35,color='w',edgecolor='k',lw=0.5)
        ax.errorbar(idx-0.2+1,mm[0],\
                    np.std(key_count_incong[key])/np.sqrt(len(key_count_incong[key])),
                    color='none',ecolor='grey',elinewidth=0.5,capsize=1,lw=0.5)
        ax.errorbar(idx+0.2+1,mm[1],\
                    np.std(key_count_congru[key])/np.sqrt(len(key_count_congru[key])),
                    color='none',ecolor='grey',elinewidth=0.5,capsize=1,lw=0.5)        
    ax.set_xticks((5,10,15,20))

    fh.savefig('4su_occurance_dist.pdf',bbox_inches='tight')
    
          
    ####################
    

def perSessSuPatt():
    
    key_count_congru=[]
    key_count_incong=[]
    
  
    
    for sess_id in range(1,114):
        
        sess_congru=0
        sess_incong=0
        if not os.path.isfile(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id)):
            print(f'missing sessid {sess_id}')
            continue
    
        filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
        if not os.path.isfile(filt_fpath):
            continue
    
        mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    
        ## count candidiates
    
        pref=np.array(mat['prefs'])[:,0]
        rego=mat['regs']
        regs=[x[0][0] for x in rego]
    
        reg_set=np.unique(regs)
        reg_set=reg_set[reg_set!='Unlabeled']
    
    
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
                sess_congru+=1
    
            else:
                sess_incong+=1
        

        key_count_congru.append(sess_congru)
        key_count_incong.append(sess_incong)
    
    # exported from matlab
    print(stats.ranksums(key_count_congru,key_count_incong))
    breakpoint()
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    
    
    (fh, ax) = plt.subplots(1, 1, figsize=(1 / 2.54, 4 / 2.54), dpi=300)

    
    ax.scatter(np.random.random(len(key_count_incong))*0.2+0.9,
                key_count_incong,s=4,c='k',alpha=0.5,edgecolors='none')
    
    ax.scatter(np.random.random(len(key_count_congru))*0.2+2.9,
               key_count_congru,s=4,c='k',alpha=0.5,edgecolors='none')
    ax.errorbar(4,np.mean(key_count_congru),\
                    np.std(key_count_congru)/np.sqrt(len(key_count_congru)),
                    fmt='ko',ecolor='k',elinewidth=0.5,capsize=2,ms=4,mfc='none') 
        
    ax.errorbar(0,np.mean(key_count_incong),\
                    np.std(key_count_incong)/np.sqrt(len(key_count_incong)),
                    fmt='ko',ecolor='k',elinewidth=0.5,capsize=2,ms=4,mfc='none')         
        
    
    ax.set_xticks([])
    ax.set_xlim([-1,5])
    ax.set_ylim([50,5000])
    ax.set_yscale('log')
    ax.set_yticks([100,1000])
    plt.show()
    breakpoint()
    fh.savefig('4su_occurance_count.pdf',bbox_inches='tight')    
    
    
def tempo_hz(congru_stats,incong_stats):
    congru_per_key={}
    incong_per_key={}
    keys=[0,1,2,3,11,12,13,22,23,33,111,112,113,122,123,133,222,223,233,333]
    for key in keys:
        congru_per_key[key]={'data':[],'mean':None,'ci':None,'n':None,'sem':None}
        incong_per_key[key]={'data':[],'mean':None,'ci':None,'n':None,'sem':None}
    
    patt_keys=[np.int64(np.dot(x.magnitude//4,[100,10,1])) for x in congru_stats['lag']]
    for idx,patt_key in enumerate(patt_keys):
        congru_per_key[patt_key]['data'].append(np.mean(congru_stats['pertrial'][idx][0]/6))

    patt_keys=[np.int64(np.dot(x.magnitude//4,[100,10,1])) for x in incong_stats['lag']]
    for idx,patt_key in enumerate(patt_keys):
        incong_per_key[patt_key]['data'].append(np.mean(incong_stats['pertrial'][idx][0]/6))


    
    for key in keys:
        congru_per_key[key]['mean']=np.mean(congru_per_key[key]['data'])
        congru_per_key[key]['ci']=boot.ci(congru_per_key[key]['data'], np.mean,n_samples=1000)
        congru_per_key[key]['n']=len(congru_per_key[key]['data'])
        congru_per_key[key]['sem']=np.std(congru_per_key[key]['data'],ddof=1)/np.sqrt(len(congru_per_key[key]['data']))
        incong_per_key[key]['mean']=np.mean(incong_per_key[key]['data'])
        incong_per_key[key]['ci']=boot.ci(incong_per_key[key]['data'], np.mean,n_samples=1000)
        incong_per_key[key]['n']=len(incong_per_key[key]['data'])
        incong_per_key[key]['sem']=np.std(incong_per_key[key]['data'],ddof=1)/np.sqrt(len(incong_per_key[key]['data']))

        
    (fh, ax) = plt.subplots(1, 1, figsize=(6 / 2.54, 4 / 2.54), dpi=300)
    
    for idx,key in enumerate(keys):
        mm=[incong_per_key[key]['mean'],congru_per_key[key]['mean']]
        
        ax.bar(idx-0.15+1,mm[0],width=0.3,color='k',edgecolor='k')
        ax.bar(idx+0.15+1,mm[1],width=0.3,color='w',edgecolor='k')
        ax.errorbar(idx-0.15+1,mm[0],\
                    incong_per_key[key]['sem'],\
                    color='none',ecolor='grey',capsize=1)
        ax.errorbar(idx+0.15+1,mm[1],\
                    congru_per_key[key]['sem'],\
                    color='none',ecolor='grey',capsize=1)        
    ax.set_xticks((5,10,15,20))
    
    
    
    ax.set_ylabel('Pattern freq. (Hz)')
    fh.savefig('4su_temproal_dist.pdf',bbox_inches='tight')
    
    
            
    (fh, ax) = plt.subplots(1, 1, figsize=(6 / 2.54, 4 / 2.54), dpi=300)
    
    for idx,key in enumerate(keys):
        mm=[incong_per_key[key]['mean'],congru_per_key[key]['mean']]
        
        ax.bar(idx-0.15+1,mm[0],width=0.3,color='k',edgecolor='k')
        ax.bar(idx+0.15+1,mm[1],width=0.3,color='w',edgecolor='k')
        ax.errorbar(idx-0.15+1,mm[0],\
                    incong_per_key[key]['sem'],\
                    color='none',ecolor='grey',capsize=1)
        ax.errorbar(idx+0.15+1,mm[1],\
                    congru_per_key[key]['sem'],\
                    color='none',ecolor='grey',capsize=1)        
    ax.set_xticks((5,10,15,20))
    
    
    
    ax.set_ylabel('Pattern freq. (Hz)')
    fh.savefig('4su_temproal_dist.pdf',bbox_inches='tight')
    
    congru_all=[]
    incong_all=[]
    for idx,key in enumerate(keys):
        congru_all.extend(congru_per_key[key]['data'])
        incong_all.extend(incong_per_key[key]['data'])
        
    
    (fh, ax) = plt.subplots(1, 1, figsize=(3 / 2.54, 4 / 2.54), dpi=300)
    
    ax.bar(1,np.mean(incong_all),width=0.8,color='k',edgecolor='k')
    ax.bar(2,np.mean(congru_all),width=0.8,color='w',edgecolor='k')
    ax.errorbar(1,np.mean(incong_all),\
                np.std(incong_all)/np.sqrt(len(incong_all)),\
                color='none',ecolor='grey',capsize=6)
    ax.errorbar(2,np.mean(congru_all),\
                np.std(congru_all)/np.sqrt(len(congru_all)),\
                color='none',ecolor='grey',capsize=6)
    ax.set_xticks([])
    ax.set_ylim(0,1.5)
    ax.set_ylabel('Average pattern freq. (Hz)')
    fh.savefig('4su_temproal_accu.pdf',bbox_inches='tight')
    print(stats.ranksums(congru_all,incong_all))
    
    
if __name__=='__main__':
    if 'congru_stats' not in dir():
        r=pickle.load(open('spade_stats.p','rb'))
        congru_stats=r['congru_stats']
        incong_stats=r['incongru_stats']
        
    # perPattRatio()
    perSessSuPatt()
    sys.exit()
    if False:
        tempo_hz(congru_stats,incong_stats)
    if True:
        perSessSuPatt()
        
    