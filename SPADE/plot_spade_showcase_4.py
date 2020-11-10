# -*- coding: utf-8 -*-
"""
 Created on Fri Oct 23 13:53:23 2020

@author: Libra
"""

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rcParams
import neo
import quantities as pq
import pickle
import scipy.io
from itertools import combinations


t_start=1000*pq.ms  #correction for missing t_start parameter during calling spade.concept_output
bin_size=4
sess_id=4

def patt_merge2(patt1,patt2):
    
    idx1=np.argsort(patt1['neurons'])
    idx2=np.argsort(patt2['neurons'])
    if not np.array_equal([patt1['neurons'][x] for x in idx1],
                          [patt2['neurons'][x] for x in idx2]):
        return False
    lag1=[np.int64(np.hstack(([0],patt1['lags'].magnitude))[x]) for x in idx1]
    lag2=[np.int64(np.hstack(([0],patt2['lags'].magnitude))[x]) for x in idx2]
    
    comb=combinations(range(4),2)
    for c in list(comb): 
        if (lag1[c[0]] > lag1[c[1]] and lag2[c[0]] < lag2[c[1]]) \
            or (lag1[c[0]] < lag1[c[1]] and lag2[c[0]] > lag2[c[1]]):
                return False
            
    return True
            
    
def patt_merge(patt):
    merged={}
    for p in patt:
        key=tuple(np.sort(p['neurons']))
        if key in merged:
            pass
        else:
            merged[key]=p
    return merged

def patt_merge3(patt):
    merged={}
    for p in patt:
        key=tuple(np.sort(p['neurons']))
        if (key not in merged) or len(p['times'])>len(merged[key]['times']):
            merged[key]=p
    return merged


def patt_merge4(patt,trialInfo,sess_id,denovo=False):
    if denovo or not os.path.isfile(f's{sess_id}cnt.p'):
        per_trial={'S1':{},'S2':{}}
        for t in range(trialInfo.shape[0]):
            print(f'merge {t}')
            win_onset=t*7000
            win_offset=t*7000+6000
            tkey='S1' if trialInfo[t,4]==4 else 'S2'
            merged={}
            for i,p in enumerate(patt):
                key=tuple(np.sort(p['neurons']))
                cnt=np.sum([x>=win_onset and x<win_offset for x in p['times'].magnitude])
                if (key not in merged) or cnt>merged[key]['cnt']:
                    merged[key]={'cnt':cnt,'patt_idx':i}
                    
            per_trial[tkey][t]=merged
        s1t=[]
        for t in per_trial['S1'].keys():
            cnt=[x['cnt'] for x in per_trial['S1'][t].values()]
            s1t.append([t,np.sum(np.array(cnt)>0),np.sum(cnt)])
    
        s2t=[]
        for t in per_trial['S2'].keys():
            cnt=[x['cnt'] for x in per_trial['S2'][t].values()]
            s2t.append([t,np.sum(np.array(cnt)>0),np.sum(cnt)])
            
        pickle.dump({'per_trial':per_trial,'s1t':s1t,'s2t':s2t},open(f's{sess_id}cnt.p','wb'))
    else:
        r=pickle.load(open(f's{sess_id}cnt.p','rb'))
        per_trial=r['per_trial']
        s1t=r['s1t']
        s2t=r['s2t']
    
    s1t=np.array(s1t)
    s1maxsu=np.max(s1t,axis=0)[1]
    s1max_freq=np.max(s1t[s1t[:,1]==s1maxsu,2])
    t1=s1t[np.logical_and(s1t[:,1]==s1maxsu,s1t[:,2]==s1max_freq),0][0]
    
    s2t=np.array(s2t)
    s2maxsu=np.max(s2t,axis=0)[1]
    s2max_freq=np.max(s2t[s2t[:,1]==s2maxsu,2])
    t2=s2t[np.logical_and(s2t[:,1]==s2maxsu,s2t[:,2]==s2max_freq),0][0]
    
    s1sel={'t1':t1,'patt_idx':[x['patt_idx'] for x in per_trial['S1'][t1].values()]}
    s2sel={'t2':t2,'patt_idx':[x['patt_idx'] for x in per_trial['S2'][t2].values()]}
    # [x['patt_idx'] for x in per_trial['S1'][t1].values()]
    return (s1sel,s2sel)



def preprocess(congru_stats,incongru_stats):
    per_sess={}
    for stats in (congru_stats,incongru_stats):
        for idx,sess in enumerate(stats['sess_ids']):
            if sess in per_sess:
                per_sess[sess].extend(stats['neu'][idx])
            else:
                per_sess[sess]=stats['neu'][idx]
            
            
    for sess in per_sess:
        per_sess[sess]=len(set(per_sess[sess]))
        
    return per_sess
            


if __name__=='__main__':
    
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['axes.linewidth'] = 0.5
    
    
    if 'congru_stats' not in dir():
        #from spade_stats.py
        fstr=pickle.load(open('spade_stats.p','rb'))
        congru_stats=fstr['congru_stats']
        incong_stats=fstr['incongru_stats']
    if False:
        per_sess_cnt=preprocess(congru_stats,incong_stats)
        sys.exit()
    # thres=0.05
    # diffsel=np.logical_and(\
    #             np.array(congru_stats['motif_pvalues'])[:,0]<thres,\
    #             np.array(congru_stats['motif_pvalues'])[:,1]<thres)
    

    filt_fpath=r'K:\code\SPADE\results\{}\winlen{}\filtered_patterns.npy'.format(sess_id,bin_size)
    r=np.load(filt_fpath,allow_pickle=True)
    patterns_raw=r[0]
    params=r[3]
    
    mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    trialInfo=mat['trialInfo']
    regs=mat['regs']
    patterns_lbl=[]
    for patt in patterns_raw:
        reg_grp=[regs[x] for x in patt['neurons']]
        if np.isin('Unlabeled',reg_grp):
            continue
        patterns_lbl.append(patt)
    
    (s1sel,s2sel)=patt_merge4(patterns_lbl,trialInfo,sess_id)

    patt_id_all=np.unique(np.hstack((s1sel['patt_idx'],s2sel['patt_idx'])))
    patt_others=np.nonzero([x not in patt_id_all for x in range(len(patterns_raw))])[0]
       
    r=neo.io.NeoMatlabIO(filename='K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    bl=r.read_block()
    spkt=bl.segments[0].spiketrains
    
    
    l_trial=s1sel['t1']*7000
    r_trial=s2sel['t2']*7000

    window_start=np.array([l_trial,r_trial])+1000

### TODO: tag patt perf    
    tag=[]
    for pidx, patt in enumerate(patterns_raw):
        ts=np.array(patt['times'].magnitude)+t_start.magnitude
        s1cnt=np.sum(np.logical_and(ts>=window_start[0],ts<window_start[0]+6000))
        s2cnt=np.sum(np.logical_and(ts>=window_start[1],ts<window_start[1]+6000))
        if s1cnt==s2cnt:
            tag.append(0)
        else:
            tag.append(1 if s1cnt>s2cnt else 2)
                
    colors=['grey','r','b']
    (fh, ax) = plt.subplots(2, 1, figsize=(20 / 2.54, 10 / 2.54), dpi=300)
    (fh2,ax2)= plt.subplots(patt_id_all.shape[0], 2, figsize=(20 / 2.54, 29 / 2.54), dpi=300)
    
    su_ids=np.sort(np.unique(np.hstack(
            ([patterns_raw[x]['neurons'] for x in patt_id_all]))))
    
    regax=[regs[x][0][0] for x in su_ids]
    
    regidx=np.argsort(regax)
    regax=[regax[x] for x in regidx]
    su_ids=su_ids[regidx]
    
    
    stp_pct=[]
    for nidx,one_id in enumerate(su_ids):
        one_spkt=spkt[one_id]
        one_stp_pct=[]
        for subidx, win_start in enumerate(window_start):
            curr_window=[win_start,win_start+6000]
            spkts=one_spkt.times.magnitude
            spk_win=spkts[np.logical_and(spkts>=curr_window[0],spkts<curr_window[1])]
            spk_stp_sel=np.zeros_like(spk_win)
            for pidx,patt in enumerate(patterns_raw):
                if one_id not in patt['neurons']:
                    continue
                in_p_idx=np.nonzero(np.array(patt['neurons'])==one_id)[0][0]
                pts=patt['times']+t_start
                pts_win=pts[np.logical_and(pts>=curr_window[0],pts<curr_window[1])].magnitude
                pts_win+=np.hstack(([0],patt['lags'].magnitude))[in_p_idx]
                patt_ts=np.full_like(spk_win,False,dtype=bool)    
                for win in pts_win:
                    pattsel=np.logical_and(spk_win>=win,spk_win<win+bin_size)
                    spk_stp_sel[pattsel]=tag[pidx]
                    patt_ts[pattsel]=True
                ##### updated f5 plot sequence
                ax[subidx].plot(\
                    np.vstack((spk_win[patt_ts],spk_win[patt_ts])),\
                    [[nidx-0.4]*np.sum(patt_ts),[nidx+0.4]*np.sum(patt_ts)],\
                    '-',color=colors[tag[pidx]],lw=0.5,alpha=1)

            ax[subidx].plot(np.vstack((spk_win[spk_stp_sel==0],spk_win[spk_stp_sel==0])),
                   [[nidx-0.4]*np.sum(spk_stp_sel==0),[nidx+0.4]*np.sum(spk_stp_sel==0)],
                       '-',color='silver',lw=0.5,alpha=1)
            one_stp_pct.append([np.sum(spk_stp_sel==0),np.sum(spk_stp_sel==1),np.sum(spk_stp_sel==2)])
                ##########################
        stp_pct.append({'yidx':nidx,'neu_id':one_id,'total_unlabeld':one_stp_pct})
    pickle.dump(stp_pct,open('stp_pct.p','wb'))
    for iidx,pidx in enumerate(patt_id_all):
        patt=patterns_raw[pidx]
        print(iidx)
        neu=patt['neurons']
        times=patt['times']
        ax2[iidx][0].set_ylabel(f'#{pidx}',fontdict={'fontsize':5})
        
        for idx,one_neu in enumerate(neu):
            for subidx, win_start in enumerate(window_start):
                curr_window=[win_start,win_start+6000]
                spkts=spkt[one_neu].times.magnitude
                spk_win=spkts[np.logical_and(spkts>=curr_window[0],spkts<curr_window[1])]
                pts=patt['times']+t_start
                pts_win=pts[np.logical_and(pts>=curr_window[0],pts<curr_window[1])].magnitude
                pts_win+=np.hstack(([0],patt['lags'].magnitude))[idx]
                
                in_stp_sel=np.full_like(spk_win,False,dtype=bool)
                
                for win in pts_win:
                    in_stp_sel[np.logical_and(spk_win>=win,spk_win<win+bin_size)]=True
                
                ax2[iidx][subidx].plot([spk_win[np.logical_not(in_stp_sel)],\
                                        spk_win[np.logical_not(in_stp_sel)]],\
                            [[idx-0.4]*np.sum(np.logical_not(in_stp_sel)),\
                             [idx+0.4]*np.sum(np.logical_not(in_stp_sel))],\
                                '-',color='silver',lw=0.25,alpha=0.5)
                if np.sum(in_stp_sel)>0:
                    ax2[iidx][subidx].plot([spk_win[in_stp_sel],spk_win[in_stp_sel]],
                                [[idx-0.4]*np.sum(in_stp_sel),[idx+0.4]*np.sum(in_stp_sel)],\
                                    '-',color=colors[tag[pidx]],lw=0.25,alpha=1)
    
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
    ax[1].set_xlabel('Tims (s), S1 trial # {}, S2 trial # {}'
                     .format(s1sel['t1']+1,s2sel['t2']+1))
    ax[1].set_yticks(range(len(su_ids)))
    ax[1].set_yticklabels(regax,fontdict={'fontsize':5})
    ax[1].set_ylim(-0.6,len(su_ids)-0.4)
    fh.savefig(f'multiple_stp_{sess_id}.pdf',bbox_inches='tight')

    for twoax in ax2:
        for idx, oneax in enumerate(twoax):
            oneax.set_xticks([])
            oneax.set_yticks([])
            oneax.set_xlim((window_start[idx],window_start[idx]+6000))
    fh2.savefig(f'expanded_stp_{sess_id}.pdf',bbox_inches='tight')
    
    (fh3,ax3) = plt.subplots(2, 1, figsize=(5 / 2.54, 10 / 2.54), dpi=300)
    for subidx,subax in enumerate(ax3):
        for neu in stp_pct:
            totalspk=np.sum(neu['total_unlabeld'][subidx])
            ratio1=neu['total_unlabeld'][subidx][1]/totalspk
            ratio2=neu['total_unlabeld'][subidx][2]/totalspk
            subax.barh(neu['yidx'],ratio1*100,color='r',edgecolor='none')
            subax.barh(neu['yidx'],ratio2*100,color='b',edgecolor='none',left=ratio1*100)
            
        subax.set_xlim(0,100)
        subax.set_yticks(range(len(su_ids)))
        subax.set_yticklabels(regax,fontdict={'fontsize':5})
        subax.set_ylim(-0.6,len(su_ids)-0.4)
    subax.set_xlabel('Spikes associated with STPs (%)')
    fh3.savefig(f'spike_ratio_{sess_id}.pdf',bbox_inches='tight')