# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 16:01:16 2020

@author: Libra
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pickle
import csv
from itertools import combinations


if __name__ == '__main__':

    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['axes.linewidth'] = 0.5

    if 'congru_stats' not in dir():
        # from spade_stats.py
        fstr = pickle.load(open('spade_stats.p', 'rb'))
        congru_stats = fstr['congru_stats']
        incong_stats = fstr['incongru_stats']

    all_reg = []
    for regs in congru_stats['neu_regs']:
        all_reg.extend(set(regs))
    all_reg = set(all_reg)

    comb_all = combinations(all_reg, 2)
    key_all = []
    for one_comb in comb_all:
        key_all.append(tuple(sorted(one_comb)))

    conn_map = {}

    for regs in congru_stats['neu_regs']:
        ureg = set(regs)
        if len(ureg) > 1:
            comb = combinations(ureg, 2)
            for c in list(comb):
                key = tuple(sorted(c))
                if key in conn_map:
                    conn_map[key] += 1
                else:
                    conn_map[key] = 1

    conn_map_fc = {}
    with open(r'K:\code\jpsth\selec_sum_ratio.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            key = tuple(sorted(row[:2]))
            conn_map_fc[key] = 1

    reg_coord = {}
    with open(r'K:\code\SPADE\reg_coord.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            key = row[0]
            reg_coord[key] = [float(row[1]), float(row[2])]

    all_reg = []
    with open('STP4_conn.csv', 'w') as f:
        write = csv.writer(f)
        write.writerow(['Source', 'Target', 'Weight', 'Partition'])
        for key in key_all:
            all_reg.extend(key)
            if key in conn_map.keys():
                # breakpoint()
                row = []
                row.extend(list(key))
                row.append(conn_map[key])
                row.append(1)
                write.writerow(row)
            elif key in conn_map_fc.keys():
                row = []
                row.extend(list(key))
                row.append(1)
                row.append(0)
                write.writerow(row)

    with open('STP4_node_coord.csv', 'w') as f:
        write = csv.writer(f)
        write.writerow(['Id', 'Label', 'AP', 'DV'])
        for reg in set(all_reg):
            row = [reg, reg, reg_coord[reg][0]/15, 1000/15-reg_coord[reg][1]/15]
            write.writerow(row)
            
        
    breakpoint()
    
    
    
def reg_degree():    
    full_reg=[]
    for regs in conn_map_fc.keys():
        full_reg.extend(set(regs))
    full_reg=set(full_reg)
    
    
    per_reg_degree={}
    for one_reg in full_reg:
        fc_cnt=0
        stp_cnt=0
        for key in conn_map_fc.keys():
            if one_reg in key:
                fc_cnt+=1
        for key in conn_map.keys():
            if one_reg in key:
                stp_cnt+=1
        per_reg_degree[one_reg]=[fc_cnt,stp_cnt]
    
        
    hfc,efc=np.histogram([x[0] for x in per_reg_degree.values()],bins=np.arange(-0.5,51))
    hstp,estp=np.histogram([x[1] for x in per_reg_degree.values()],bins=np.arange(-0.5,51))

    (fh,ax) = plt.subplots(1, 1, figsize=(5 / 2.54, 5 / 2.54), dpi=300)
    ax.plot(efc[1:]-0.5,hfc,'-',color='silver')
    ax.plot(estp[1:]-0.5,hstp,'-r')
    ax.set_ylim(0,20)
    ax.set_xlabel('Degree')
    ax.set_ylabel('Number of region')
    fh.save('reg_degree.pdf',bbox_inches='tight')
    
                
            
