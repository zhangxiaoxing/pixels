# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 13:31:44 2020

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
# import matlab.engine
import csv
import h5py

denovo=True
t_start=1000*pq.ms  #correction for missing t_start parameter during calling spade.concept_output
bin_size=4

stp_pct=[]
if __name__=='__main__':
    # meng=matlab.engine.start_matlab('-nodesktop -sd "K:\code\jpsth"')
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['axes.linewidth'] = 0.5
    
        
    if denovo:
        spk_ts_export=[]
        spk_count=[]
        for sess_id in range(114):
            print(f'sess {sess_id}')
            spkf=r'K:\code\SPADE\spkt\spktN13_{}.mat'.format(sess_id)
            if not os.path.isfile(spkf):
                continue
            mat=scipy.io.loadmat(spkf)
            regs=mat['regs']
            samp=mat['trialInfo'][:,4]
            
            r=neo.io.NeoMatlabIO(filename=spkf.replace('N13', 'O17'))
            bl=r.read_block()
            spkt=bl.segments[0].spiketrains
    
            filt_fpath=r'K:\code\SPADE\results\{}\winlen{}\filtered_patterns.npy'.format(sess_id,bin_size)
            if not os.path.isfile(filt_fpath):
                reg_sel=[x[0][0]!='Unlabeled' for x in regs]
                spkt=[spkt[x] for x in reg_sel if x]
                stp_pct.extend(np.zeros(len(spkt)))
                continue
            
            
            r=np.load(filt_fpath,allow_pickle=True)
            patterns_raw=r[0]
            params=r[3]
            
            patt_neu_set=[]
            for patt in patterns_raw:
                reg_grp=[regs[x] for x in patt['neurons']]
                if np.isin('Unlabeled',reg_grp):
                    continue
                patt_neu_set.extend(patt['neurons'])
    
            
            patt_neu_all=list(set(patt_neu_set))
            patt_others=np.nonzero([x not in patt_neu_all for x in range(len(spkt))])[0]
            
            stp_pct.extend(np.zeros_like(patt_others))
               
            for nidx,one_id in enumerate(patt_neu_all):
                one_spkt=spkt[one_id]
                one_stp_pct=[]
                spkts=one_spkt.times.magnitude
                spk_stp_sel=np.zeros_like(spkts)
                neu_win=[]
                for pidx,patt in enumerate(patterns_raw):
                    if one_id not in patt['neurons']:
                        continue
                    in_p_idx=np.nonzero(np.array(patt['neurons'])==one_id)[0][0]
                    pts=patt['times'].magnitude+t_start.magnitude
                    pts+=np.hstack(([0],patt['lags'].magnitude))[in_p_idx]
                    neu_win.extend(pts)
                for win in list(set(neu_win)):
                    pattsel=np.logical_and(spkts>=win,spkts<win+bin_size)
                    spk_stp_sel[pattsel]=1
### Export for FC tag
                if np.sum(spk_stp_sel>0)>0:
                    spt_ts=spkts[spk_stp_sel>0]
                    for ts in spt_ts:
                        row = [sess_id,mat['transIds'][nidx][0],samp[int(ts//7000)],ts]
                        spk_ts_export.append(row)
    
                stp_pct.append(np.sum(spk_stp_sel)/spk_stp_sel.shape[0])
                spk_count.append([sess_id,mat['transIds'][nidx][0],np.sum(spk_stp_sel),spk_stp_sel.shape[0]])
        
        with h5py.File('stp_ts_export.hdf5','w') as fw:
            fw.create_dataset('spk_ts',data=np.array(spk_ts_export))
            fw.create_dataset('spk_cnt',data=np.array(spk_count))
        pickle.dump(stp_pct,open('stp_pct_all.p','wb'))
    else:
        stp_pct=pickle.load(open('stp_pct_all.p','rb'))
        
    breakpoint()
        
    count=[]
    for e in np.arange(-0.1,1,0.1):
        count.append(np.sum([x>e and x<=e+0.1 for x in stp_pct]))
    
    # meng.quit()
    (fh, ax) = plt.subplots(1, 1, figsize=(5 / 2.54, 5 / 2.54), dpi=300)
    ax.bar(np.arange(0.05,1,0.1)*100,count[1:]/np.sum(count[1:]),width=8,color='w',edgecolor='k')
    ax.set_ylabel('Number of neurons')
    ax.set_xlabel('STP associated spikes (%)')
    # ax.set_ylim(0,40)
    fh.savefig('STP_asso_spks.pdf',bbox_inches='tight')
    #total memory neuron is len(stp_pct)
        
