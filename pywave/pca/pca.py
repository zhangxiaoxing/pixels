# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:57:04 2021

@author: Libra
"""

import numpy as np
import pywave.datautils as util
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, FastICA
from sklearn.preprocessing import StandardScaler


def pca_3d(delay=6):
    fr_per_su = util.get_dataset()
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])

    scaler=StandardScaler()
    delaybin=range(16,16+delay*4);
    scaler.fit(np.vstack((S1_mm[delaybin,:],S2_mm[delaybin,:])))
    norm_mm=scaler.transform(np.vstack((S1_mm,S2_mm)))

    pca = PCA(n_components=20)
    #TODO 3s delay compatbility
    pcamat=norm_mm[list(delaybin)+list(np.array(delaybin)+56),:]
    pca.fit(pcamat)
    comp=pca.transform(norm_mm)
    coeff = pca.components_
    ratio = pca.explained_variance_ratio_

    fig=plt.figure()
    ax = fig.add_subplot(projection='3d')
    h1=ax.plot(comp[8:45,0],comp[8:45,1],comp[8:45,2],'-r')
    h2=ax.plot(comp[64:101,0],comp[64:101,1],comp[64:101,2],'-b')

    interp=[np.interp(np.array([12,16,12+56,16+56,40,44,40+56,44+56])-0.5,
                      range(comp.shape[0]),
                      comp[:,c])
            for c in range(3)]
    hs=[ax.plot(interp[0][t],interp[1][t],interp[2][t],'o',c=[0.5,0.5,0.5]) for t in range(4)]
    ht=[ax.plot(interp[0][t],interp[1][t],interp[2][t],'ko') for t in range(4,8)]

    ax.set_zlabel('PC3')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    plt.legend([h1[0],h2[0],hs[0][0],ht[0][0]], ['S1 trial','S2 trial','Sample','Test'])


def pca_dist(delay=6):
    fr_per_su = util.get_dataset()
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])

    scaler=StandardScaler()
    delaybin=range(16,16+delay*4);
    scaler.fit(np.vstack((S1_mm[delaybin,:],S2_mm[delaybin,:])))
    norm_mm=scaler.transform(np.vstack((S1_mm,S2_mm)))

    pca = PCA(n_components=20)
    #TODO 3s delay compatbility
    pcamat=norm_mm[list(delaybin)+list(np.array(delaybin)+56),:]
    pca.fit(pcamat)
    comp=pca.transform(norm_mm)
    coeff = pca.components_
    ratio = pca.explained_variance_ratio_

    pcnum=8
    hp=[None]*pcnum
    plt.figure()
    for c in range(pcnum):
        hp[c]=plt.plot(np.abs(np.array(comp[:56,c])-np.array(comp[56:,c])))

#TODO: eculidian distance
    # eculidian=[np.linalg.norm(
    #     [proj1E[t]-proj2E[t],proj1M[t]-proj2M[t],proj1L[t]-proj2L[t]]) for t in range(len(proj1E))]


    plt.legend([hp[c][0] for c in range(pcnum)], [f'PC{c+1}, {ratio[c]*100:.1f}%' for c in range(pcnum)])
    ax=plt.gca()
    ax.set_ylabel('S1-S2 PC space distance (a.u.)')
    ax.set_xticks(np.array([12,32,53])-0.5)
    ax.set_xticklabels([0,5,10])
    ax.set_xlabel('Time (s)')
    testmark=[40,44] if delay==6 else [28,32]
    [plt.axvline(x,ls='--',c='k') for x in np.array([12,16]+testmark)-0.5]
    plt.title(f'PC1-PC{pcnum}, total {np.sum(ratio[:pcnum])*100:.1f}% delay variance')

def ica_dist(delay=6):
    fr_per_su = util.get_dataset()
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])

    scaler=StandardScaler()
    delaybin=range(16,16+delay*4);
    scaler.fit(np.vstack((S1_mm[delaybin,:],S2_mm[delaybin,:])))
    norm_mm=scaler.transform(np.vstack((S1_mm,S2_mm)))

    ica = FastICA(n_components=20)
    #TODO 3s delay compatbility
    pcamat=norm_mm[list(delaybin)+list(np.array(delaybin)+56),:]
    ica.fit(pcamat)
    comp=ica.transform(norm_mm)
    coeff = ica.mixing_

    icnum=5
    hp=[None]*icnum
    plt.figure()
    for c in range(icnum):
        hp[c]=plt.plot(np.abs(np.array(comp[:56,c])-np.array(comp[56:,c])))

    # eculidian=[np.linalg.norm([proj1E[t]-proj2E[t],proj1M[t]-proj2M[t],proj1L[t]-proj2L[t]]) for t in range(len(proj1E))]
#
    plt.legend([hp[c][0] for c in range(icnum)], [f'IC{c+1}' for c in range(icnum)])
    ax=plt.gca()
    ax.set_ylabel('S1-S2 PC space distance (a.u.)')
    ax.set_xticks(np.array([12,32,53])-0.5)
    ax.set_xticklabels([0,5,10])
    ax.set_xlabel('Time (s)')
    testmark=[40,44] if delay==6 else [28,32]
    [plt.axvline(x,ls='--',c='k') for x in np.array([12,16]+testmark)-0.5]
    plt.title(f'IC1-IC{icnum}')


    fig=plt.figure()
    ax = fig.add_subplot(projection='3d')
    h1=ax.plot(comp[8:45,0],comp[8:45,1],comp[8:45,2],'-r')
    h2=ax.plot(comp[64:101,0],comp[64:101,1],comp[64:101,2],'-b')

    interp=[np.interp(np.array([12,16,12+56,16+56,40,44,40+56,44+56])-0.5,
                      range(comp.shape[0]),
                      comp[:,c])
            for c in range(3)]
    hs=[ax.plot(interp[0][t],interp[1][t],interp[2][t],'o',c=[0.5,0.5,0.5]) for t in range(4)]
    ht=[ax.plot(interp[0][t],interp[1][t],interp[2][t],'ko') for t in range(4,8)]

    ax.set_zlabel('IC3')
    ax.set_xlabel('IC1')
    ax.set_ylabel('IC2')
    plt.legend([h1[0],h2[0],hs[0][0],ht[0][0]], ['S1 trial','S2 trial','Sample','Test'])