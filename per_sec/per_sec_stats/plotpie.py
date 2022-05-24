# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:58:12 2021

@author: Libra
"""
import matplotlib.pyplot as plt
from matplotlib import rcParams


def plotpie():
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    # frac = [switched_count + unclassified_count + bs_count + sample_only_count, sust_count, transient_count,
    #         non_sel_mod_count + non_mod_count, ]
    frac = [sust_count, transient_count+unclassified_count,
            non_sel_mod_count + non_mod_count + switched_count + bs_count + sample_only_count,]
    print(np.sum(frac))
    explode = (0.1, 0.1, 0)
    labels = ('sustained', 'transient', 'non-selective')

    (fh, ax) = plt.subplots(1, 1, figsize=(12 / 2.54, 4 / 2.54), dpi=300)
    if counterclock:
        startangle = -60
    else:
        startangle = 180
    ax.pie(frac, explode=explode, labels=labels, autopct='%1.1f%%',
           shadow=False, startangle=startangle,counterclock=counterclock,
           colors=('blue', 'red', 'w'),
           wedgeprops={'ls':'-','lw':0.5,'ec':'k'})
    # ax.axis('equal')
    ax.set_xlim((-0.7, 0.7))
    # fh.suptitle(f'{delay}s delay')
    plt.show()
    fh.savefig(f'sus_trans_pie_{delay}.pdf')