# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 13:53:23 2020

@author: Libra
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import neo
import quantities as pq
import pickle
import scipy.io


t_start=1000*pq.ms  #correction for missing t_start parameter during calling spade.concept_output
bin_size=4

s1idx=5
s2idx=5

if __name__=='__main__':
    
    if 'congru_stats' not in dir():
        #from spade_stats.py
        fstr=pickle.load(open('spade_stats.p','rb'))
        congru_stats=fstr['congru_stats']
        incong_stats=fstr['incongru_stats']
    
    thres=0.05
    diffsel=np.logical_and(\
                np.array(congru_stats['motif_pvalues'])[:,0]<thres,\
                np.array(congru_stats['motif_pvalues'])[:,1]<thres)
    
    # sess_id=np.argmax(np.histogram(
    #     np.array(congru_stats['sess_ids'])[diffsel],range(114))[0])
    sess_id=7

    filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
    r=np.load(filt_fpath,allow_pickle=True)
    patterns=r[0]
    params=r[3]

    # plot_idces=np.logical_and(diffsel,np.array(congru_stats['sess_ids'])==sess_id)
    plot_idces=np.array(congru_stats['sess_ids'])==sess_id
        
    r=neo.io.NeoMatlabIO(filename='K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    bl=r.read_block()
    spkt=bl.segments[0].spiketrains
    
    mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    trialInfo=mat['trialInfo']
    s1_trial=np.nonzero(trialInfo[:,4]==4)[0]
    s2_trial=np.nonzero(trialInfo[:,4]==8)[0]
    regs=mat['regs']
    
    l_trial=s1_trial[s1idx]*7000
    r_trial=s2_trial[s2idx]*7000

    window_start=np.array([l_trial,r_trial])+1000

    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=np.sum(plot_idces)-1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


    (fh, ax) = plt.subplots(2, 1, figsize=(20 / 2.54, 10 / 2.54), dpi=300)
    
#--- TODO: plot gray background    
    for nidx,one_spkt in enumerate(spkt):
        for subidx, win_start in enumerate(window_start):
            curr_window=[win_start,win_start+6000]
            spkts=one_spkt.times.magnitude
            spk_win=spkts[np.logical_and(spkts>=curr_window[0],spkts<curr_window[1])]

            ax[subidx].plot(np.vstack((spk_win,spk_win)),\
                   [[nidx-0.5]*spk_win.shape[0],[nidx+0.5]*spk_win.shape[0]],\
                       '-',color='silver',lw=0.5,alpha=0.5)


    
    
    for iidx,pre_conversion_idx in enumerate(np.nonzero(plot_idces)[0]):
        neu=congru_stats['neu'][pre_conversion_idx]
        sib=np.nonzero([np.array_equal(neu,p['neurons']) for p in patterns])
        sigsib=np.nonzero([(congru_stats['sess_ids'][i]==sess_id 
                            and np.array_equal(neu,congru_stats['neu'][i]))
                            for i in range(len(congru_stats['sess_ids']))])[0]
        patt_plot_idx=sib[0][np.nonzero(sigsib==pre_conversion_idx)[0][0]]
        
        times=patterns[patt_plot_idx]['times']
        
        colorVal = scalarMap.to_rgba(iidx)
        
        for idx,one_neu in enumerate(patterns[patt_plot_idx]['neurons']):
            for subidx, win_start in enumerate(window_start):
                curr_window=[win_start,win_start+6000]
                spkts=spkt[one_neu].times.magnitude
                spk_win=spkts[np.logical_and(spkts>=curr_window[0],spkts<curr_window[1])]
                
                pts=patterns[patt_plot_idx]['times']+t_start
                pts_win=pts[np.logical_and(pts>=curr_window[0],pts<curr_window[1])].magnitude
                pts_win+=np.hstack(([0],patterns[patt_plot_idx]['lags'].magnitude))[idx]
                
                in_stp_sel=np.full_like(spk_win,False,dtype=bool)
                
                for win in pts_win:
                    in_stp_sel[np.logical_and(spk_win>=win,spk_win<win+bin_size)]=True
                
                # ax.plot(spk_win[np.logical_not(in_stp_sel)],\
                        # [one_neu]*np.sum(np.logical_not(in_stp_sel)),'|',color='silver')
                if np.sum(in_stp_sel)>0:
                    ax[subidx].plot([spk_win[in_stp_sel],spk_win[in_stp_sel]],
                                [[one_neu-0.5]*np.sum(in_stp_sel),[one_neu+0.5]*np.sum(in_stp_sel)],\
                                    '-',color=colorVal,lw=1,alpha=1)
        
    
    ax[0].set_xlim((window_start[0],window_start[0]+6000))
    ax[0].set_yticks(range(0,len(spkt),5))
    ax[0].set_yticklabels([regs[x][0][0] for x in range(0,len(spkt),5)])
    ax[0].set_xticks(np.arange(window_start[0],window_start[0]+6001,1000))
    ax[0].set_xticklabels(np.arange(7))
    ax[0].set_xlabel('Tims (s), S1 trial # {}'.format(s1_trial[s1idx]+1))
    ax[1].set_xlim((window_start[1],window_start[1]+6000))
    ax[1].set_xticks(np.arange(window_start[1],window_start[1]+6001,1000))
    ax[1].set_xticklabels(np.arange(7))
    ax[1].set_xlabel('Tims (s), S2 trial # {}'.format(s2_trial[s2idx]+1))
