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
    by_neu={}
    if 'congru_stats' not in dir():
        r=pickle.load(open('spade_stats.p','rb'))
        congru_stats=r['congru_stats']
    if True:
        for idx,neu in enumerate(congru_stats['neu']):
            # breakpoint()
            sneu=tuple(sorted(neu))
            if sneu in by_neu:
                by_neu[sneu].append([idx,neu,congru_stats['lag'][idx],congru_stats['pertrial'][idx]])
            else:
                by_neu[sneu]=[[idx,neu,congru_stats['lag'][idx],congru_stats['pertrial'][idx]]]
                       
        all_sneu=[]
        for sneu in by_neu.values():
            maxs1=[0,0,0,0,0] #idx, patt_id, neu, trial, count
            maxs2=[0,0,0,0,0] #idx, patt_id, neu, trial, count
            
            for idx,patt in enumerate(sneu):
                if not np.array_equal(patt[2],[4,8,12]):
                    continue
                m1=np.max(patt[3][0])
                if m1> maxs1[4]:
                    maxs1=[idx,patt[0],patt[1],np.argmax(patt[3][0]),m1]
                m2=np.max(patt[3][1])
                if m2> maxs2[4]:
                    maxs2=[idx,patt[0],patt[1],np.argmax(patt[3][1]),m2]
            
            all_sneu.append([maxs1,maxs2])
                    
            
        # for idx,sneu in enumerate(all_sneu):
        #     if not np.array_equal(sneu[0][1],sneu[1][1]):
        #         print([idx,sneu[0][1],sneu[1][1]])
    if False:
        reg_of_4=np.nonzero([np.unique(x).shape[0]==4 for x in congru_stats['neu_regs']])[0]
        for i in reg_of_4:
            print(congru_stats['neu_regs'][i]) #3447

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
    
    # for sneu_idx in [33,47,77,78,84]:        
    #     patt_idx_all=[all_sneu[sneu_idx][0][1],all_sneu[sneu_idx][1][1]]
        
    #     sess_id=congru_stats['sess_ids'][patt_idx_all[0]]
    #     l_idx=all_sneu[sneu_idx][0][3]
    #     r_idx=all_sneu[sneu_idx][1][3]
    
    for pidx in [3447]:
        sess_id=congru_stats['sess_ids'][pidx]

        
        mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
        trialInfo=mat['trialInfo']
        s1sel=trialInfo[:,4]==4
        s2sel=trialInfo[:,4]==8
        
        goodsel=np.logical_and(trialInfo[:,8],trialInfo[:,9])
        
        l_trial=np.nonzero(np.logical_and(s1sel,goodsel))[0][l_idx]*7000
        r_trial=np.nonzero(np.logical_and(s2sel,goodsel))[0][r_idx]*7000
    
        window_start=np.array([l_trial,r_trial])+1000
    
        rcParams['pdf.fonttype'] = 42
        rcParams['ps.fonttype'] = 42
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Arial']
    
        colors=['grey','r','b']
        (fh, ax) = plt.subplots(2, 1, figsize=(20 / 2.54, 4 / 2.54), dpi=300)
        
        su_ids=congru_stats['neu'][patt_idx_all[0]]
        
        regax=congru_stats['neu_regs'][patt_idx_all[0]]
    
        filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
        r=np.load(filt_fpath,allow_pickle=True)
        patterns=r[0]
            
        r=neo.io.NeoMatlabIO(filename='K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
        bl=r.read_block()
        spkt=bl.segments[0].spiketrains
    
        
        # stp_pct=[]
        for nidx,one_id in enumerate(su_ids):
            one_spkt=spkt[one_id]
            # one_stp_pct=[]
            for subidx, win_start in enumerate(window_start):
                curr_window=[win_start,win_start+6000]
                spkts=one_spkt.times.magnitude
                spk_win=spkts[np.logical_and(spkts>=curr_window[0],spkts<curr_window[1])]
                spk_stp_sel=np.zeros_like(spk_win)
                for pidx,patt_idx in enumerate(patt_idx_all):
                    
                    neu=congru_stats['neu'][patt_idx]
                    # patt_reg=congru_stats['neu_regs'][patt_idx]
                    sib=np.nonzero([np.array_equal(neu,p['neurons']) for p in patterns])
                    sigsib=np.nonzero([(congru_stats['sess_ids'][i]==sess_id 
                                        and np.array_equal(neu,congru_stats['neu'][i]))
                                        for i in range(len(congru_stats['sess_ids']))])[0]
                    
                    patt_plot_idx=sib[0][np.nonzero(sigsib==patt_idx)[0][0]]
                    
                    in_p_idx=np.nonzero(np.array(congru_stats['neu'][patt_idx])==one_id)[0][0]
                    pts=patterns[patt_plot_idx]['times']+t_start
                    pts_win=pts[np.logical_and(pts>=curr_window[0],pts<curr_window[1])].magnitude
                    pts_win+=np.hstack(([0],patterns[patt_plot_idx]['lags'].magnitude))[in_p_idx]
                    patt_ts=np.full_like(spk_win,False,dtype=bool)    
                    for win in pts_win:
                        pattsel=np.logical_and(spk_win>=win,spk_win<win+bin_size)
                        spk_stp_sel[pattsel]=pidx+1
                        patt_ts[pattsel]=True
                    ##### updated f5 plot sequence
                    ax[subidx].plot(\
                        np.vstack((spk_win[patt_ts],spk_win[patt_ts])),\
                        [[nidx-0.4]*np.sum(patt_ts),[nidx+0.4]*np.sum(patt_ts)],\
                        '-',color=colors[pidx+1],lw=0.25,alpha=0.75)
    
                ax[subidx].plot(np.vstack((spk_win[spk_stp_sel==0],spk_win[spk_stp_sel==0])),
                       [[nidx-0.4]*np.sum(spk_stp_sel==0),[nidx+0.4]*np.sum(spk_stp_sel==0)],
                           '-',color='k',lw=0.25,alpha=0.15)
                # one_stp_pct.append([np.sum(spk_stp_sel==0),np.sum(spk_stp_sel==1),np.sum(spk_stp_sel==2)])
            # stp_pct.append({'yidx':nidx,'neu_id':one_id,'total_unlabeld':one_stp_pct})
                    ##########################
        
        
        ax[0].set_xlim((window_start[0],window_start[0]+6000))
        ax[0].set_yticks(range(len(su_ids)))
        ax[0].set_yticklabels(regax,fontdict={'fontsize':5})
        ax[0].set_xticks(np.arange(window_start[0],window_start[ 0]+6001,1000))
        ax[0].set_xticklabels(np.arange(7))
        # ax[0].set_xlabel('Tims (s), S1 trial # {}'.format(s1sel['t1']+1))
        ax[0].set_ylim(-0.6,len(su_ids)-0.4)
        ax[1].set_xlim((window_start[1],window_start[1]+6000))
        ax[1].set_xticks(np.arange(window_start[1],window_start[1]+6001,1000))
        ax[1].set_xticklabels(np.arange(7))
        ax[1].set_yticks(range(4))
        ax[1].set_yticklabels(regax,fontdict={'fontsize':5})
        ax[1].set_xlabel('Mouse #  ,Sess {}, S1 trial #{}, S2 trial #{} Tims (s)'\
                         .format(sess_id,
                                 np.nonzero(np.logical_and(s1sel,goodsel))[0][l_idx],
                                 np.nonzero(np.logical_and(s1sel,goodsel))[0][r_idx]))
        plt.subplots_adjust(hspace=1)
        fh.savefig(f'4su_1patt_showcase_{sneu_idx}.pdf',bbox_inches='tight')
        plt.close('all')
    sys.exit()





#####################Previous plot script
def previous_plot():
    (fh, ax) = plt.subplots(2, 1, figsize=(20 / 2.54, 5 / 2.54), dpi=300)

    
    for patt_idx in patt_idx_all:
    
        neus=congru_stats['neu'][patt_idx]
        
        filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
        r=np.load(filt_fpath,allow_pickle=True)
        patterns=r[0]
    
        plot_idces=patt_idx
            
        r=neo.io.NeoMatlabIO(filename='K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
        bl=r.read_block()
        spkt=bl.segments[0].spiketrains
        
        neu=congru_stats['neu'][patt_idx]
        patt_reg=congru_stats['neu_regs'][patt_idx]
        sib=np.nonzero([np.array_equal(neu,p['neurons']) for p in patterns])
        sigsib=np.nonzero([(congru_stats['sess_ids'][i]==sess_id 
                            and np.array_equal(neu,congru_stats['neu'][i]))
                            for i in range(len(congru_stats['sess_ids']))])[0]
        
        patt_plot_idx=sib[0][np.nonzero(sigsib==patt_idx)[0][0]]
        times=patterns[patt_plot_idx]['times']
        
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
        # fh.savefig('4su_1patt_showcase.pdf',bbox_inches='tight')
        