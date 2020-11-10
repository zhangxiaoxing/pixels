# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 08:45:43 2020

@author: Libra
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import neo
import quantities as pq
import pickle
import scipy.io
from matplotlib import rcParams


t_start=1000*pq.ms  #correction for missing t_start parameter during calling spade.concept_output
bin_size=4

if __name__=='__main__':
    if 'congru_stats' not in dir():
        r=pickle.load(open('spade_stats.p','rb'))
        congru_stats=r['congru_stats']
    if False:
        patt_match=[np.array_equal(x,[4,8,12]) for x in congru_stats['lag']]
        s_sig=[x<0.01 for x in list(zip(*congru_stats['motif_pvalues']))[0]]
        indices,=np.nonzero(np.logical_and(patt_match,s_sig))
        for idx in indices:
            print([idx,
                   np.max(congru_stats['pertrial'][idx][0]),
                   np.min(congru_stats['pertrial'][idx][0]),
                   np.max(congru_stats['pertrial'][idx][1]),
                   np.min(congru_stats['pertrial'][idx][1]),])
    
    patt_idx=4536
    maxfq=12
    neus=congru_stats['neu'][patt_idx]
    # l_trial_idx=np.nonzero(congru_stats['pertrial'][patt_idx][0]==12)[0][0]
    # r_trial_idx=np.nonzero(congru_stats['pertrial'][patt_idx][1]==0)[0][0]
    sess_id=congru_stats['sess_ids'][patt_idx]
    
    filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
    r=np.load(filt_fpath,allow_pickle=True)
    patterns=r[0]
    params=r[3]

    plot_idces=patt_idx
        
    r=neo.io.NeoMatlabIO(filename='K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    bl=r.read_block()
    spkt=bl.segments[0].spiketrains
    
    mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    trialInfo=mat['trialInfo']
    s1sel=trialInfo[:,4]==4
    s2sel=trialInfo[:,4]==8
    
    goodsel=np.logical_and(trialInfo[:,8],trialInfo[:,9])
    errorsel=np.logical_not(trialInfo[:,9])
    
    neu=congru_stats['neu'][patt_idx]
    patt_reg=congru_stats['neu_regs'][patt_idx]
    sib=np.nonzero([np.array_equal(neu,p['neurons']) for p in patterns])
    sigsib=np.nonzero([(congru_stats['sess_ids'][i]==sess_id 
                        and np.array_equal(neu,congru_stats['neu'][i]))
                        for i in range(len(congru_stats['sess_ids']))])[0]
    
    patt_plot_idx=sib[0][np.nonzero(sigsib==patt_idx)[0][0]]
    times=patterns[patt_plot_idx]['times']
    hh=np.histogram(times,np.arange(0,7000*trialInfo.shape[0]+1,7000))[0]
    l_idx=np.nonzero(np.logical_and(s1sel,hh==np.max(hh[s1sel])))[0][0]
    r_idx=np.nonzero(np.logical_and(s2sel,hh==np.min(hh[s2sel])))[0][0]

    l_trial=l_idx*7000
    r_trial=r_idx*7000

    window_start=np.array([l_trial,r_trial])+1000

    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    (fh, ax) = plt.subplots(2, 1, figsize=(20 / 2.54, 5 / 2.54), dpi=300)
    
#--- TODO: plot gray background    
    for nidx,one_spkt in enumerate([spkt[x] for x in neus]):
        for subidx, win_start in enumerate(window_start):
            curr_window=[win_start,win_start+6000]
            spkts=one_spkt.times.magnitude
            spk_win=spkts[np.logical_and(spkts>=curr_window[0],spkts<curr_window[1])]

            ax[subidx].plot(np.vstack((spk_win,spk_win)),\
                   [[nidx-0.5]*spk_win.shape[0],[nidx+0.5]*spk_win.shape[0]],\
                       '-',color='silver',lw=0.5,alpha=0.5)
    
    for idx,one_neu in enumerate(patterns[patt_plot_idx]['neurons']):
        for subidx, win_start in enumerate(window_start):
            curr_window=[win_start,win_start+6000]
            spkts=spkt[one_neu].times.magnitude
            spk_win=spkts[np.logical_and(spkts>=curr_window[0],spkts<curr_window[1])]
            
            pts=times+t_start
            pts_win=pts[np.logical_and(pts>=curr_window[0],pts<curr_window[1])].magnitude
            pts_win+=np.hstack(([0],patterns[patt_plot_idx]['lags'].magnitude))[idx]
            
            in_stp_sel=np.full_like(spk_win,False,dtype=bool)
            
            for win in pts_win:
                in_stp_sel[np.logical_and(spk_win>=win,spk_win<win+bin_size)]=True
            
            # ax.plot(spk_win[np.logical_not(in_stp_sel)],\
                    # [one_neu]*np.sum(np.logical_not(in_stp_sel)),'|',color='silver')
            if np.sum(in_stp_sel)>0:
                ax[subidx].plot([spk_win[in_stp_sel],spk_win[in_stp_sel]],
                            [[idx-0.5]*np.sum(in_stp_sel),[idx+0.5]*np.sum(in_stp_sel)],\
                                '-',color='red',lw=0.5,alpha=1)
        
    
    ax[0].set_xlim((window_start[0],window_start[0]+6000))
    ax[0].set_yticks(range(4))
    ax[0].set_yticklabels(patt_reg)
    ax[0].set_xticks(np.arange(window_start[0],window_start[0]+6001,1000))
    ax[0].set_xticklabels(np.arange(7))
    ax[0].set_xlabel('Tims (s)')
    ax[0].set_ylabel('S1 trial #{}'.format(l_idx+1))
    ax[1].set_xlim((window_start[1],window_start[1]+6000))
    ax[1].set_xticks(np.arange(window_start[1],window_start[1]+6001,1000))
    ax[1].set_xticklabels(np.arange(7))
    ax[1].set_yticks(range(4))
    ax[1].set_yticklabels(patt_reg)
    ax[1].set_xlabel(f'Sess {sess_id}, Tims (s)')
    ax[1].set_ylabel('S2 trial #{}'.format(r_idx+1))
    plt.subplots_adjust(hspace=1)
    fh.savefig('4su_1patt_showcase.pdf',bbox_inches='tight')
    