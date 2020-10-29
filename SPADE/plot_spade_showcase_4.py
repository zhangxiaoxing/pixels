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
from spade_stats import results

t_start=1000*pq.ms  #correction for missing t_start parameter during calling spade.concept_output
bin_size=4

if __name__=='__main__':
    
    if 'congru_stats' not in dir():
        #from spade_stats.py
        fstr=pickle.load(open('spade_stats.p','rb'))
        congru_stats=fstr['congru_stats']
    
    thres=0.005
    diffsel=np.nonzero(np.logical_and(\
                np.array(congru_stats.perHz_pvalues)[:,0]<thres,\
                np.array(congru_stats.perHz_pvalues)[:,1]<thres))
    
    idx=diffsel[0][0]
    sess_id=congru_stats.sess_ids[idx]
    neu=congru_stats.neu[idx]
    
    
    filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
    r=np.load(filt_fpath,allow_pickle=True)
    patterns=r[0]
    params=r[3]
    
    sib=np.nonzero([np.array_equal(neu,p['neurons']) for p in patterns])
    
    sigsib=np.nonzero([(congru_stats.sess_ids[i]==sess_id 
                        and np.array_equal(neu,congru_stats.neu[i]))
                       for i in range(len(congru_stats.sess_ids))])[0]
    patt_plot_idx=sib[0][np.nonzero(sigsib==idx)[0][0]]
    
    r=neo.io.NeoMatlabIO(filename='K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
    bl=r.read_block()
    spkt=bl.segments[0].spiketrains
    
    
    times=patterns[patt_plot_idx]['times']
    # hist,edge=np.histogram(times,np.arange(0,7000*121,7000))
    
    #--- plot
    (fh, ax) = plt.subplots(1, 1, figsize=(12 / 2.54, 4 / 2.54), dpi=300)
    for idx,one_neu in enumerate(patterns[patt_plot_idx]['neurons']):
        curr_window=[1000,7000]
        spkts=spkt[one_neu].times.magnitude
        spk_win=spkts[np.logical_and(spkts>=curr_window[0],spkts<curr_window[1])]
        
        pts=patterns[patt_plot_idx]['times']+t_start
        pts_win=pts[np.logical_and(pts>=curr_window[0],pts<curr_window[1])].magnitude
        pts_win+=np.hstack(([0],patterns[patt_plot_idx]['lags'].magnitude))[idx]
        
        in_stp_sel=np.full_like(spk_win,False,dtype=bool)
        
        for win in pts_win:
            in_stp_sel[np.logical_and(spk_win>=win,spk_win<win+bin_size)]=True
        
        ax.plot(spk_win[np.logical_not(in_stp_sel)],\
                [one_neu]*np.sum(np.logical_not(in_stp_sel)),'|',color='silver')
        ax.plot(spk_win[in_stp_sel],[one_neu]*np.sum(in_stp_sel),'|',color='red')
        
    
    ax.vlines(np.arange(1000,70000,7000),0,80,'k',':')
    ax.set_xlim((1000,2000))
