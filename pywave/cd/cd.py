# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:57:04 2021

@author: Libra
"""

import numpy as np
from pywave import datautils as util
import matplotlib.pyplot as plt


def cd_heatmap(delay=6,bound=0.5):
    fr_per_su = util.get_dataset()
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])

    cdMat=(S1_mm-S2_mm)
    cdco=np.corrcoef(cdMat)

    plt.figure()
    im = plt.imshow(cdco, cmap="jet", aspect="equal", vmin=-1, vmax=1,origin='lower')
    testmark=[40,44] if delay==6 else [28,32]
    [plt.axhline(x,ls='--',c='w') for x in np.array([12,16]+testmark)-0.5]
    [plt.axvline(x,ls='--',c='w') for x in np.array([12,16]+testmark)-0.5]
    plt.colorbar(im)
    plt.title('C.D. cross-temporal corr. coef.')
    ax=plt.gca()
    ax.set_xticks([12.5, 32.5, 52.5])
    ax.set_xticklabels([0, 5, 10])
    ax.set_yticks([12.5, 32.5, 52.5])
    ax.set_yticklabels([0, 5, 10])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Time (s)')
    ax.set_xlim(7.5,43.5)
    ax.set_ylim(7.5,43.5)



def cdT_projection(delay=6):
    fr_per_su = util.get_dataset()
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])

    cdMat=(S1_mm-S2_mm)
    cdDelayE=np.mean(cdMat[16:24,:],axis=0)
    cdDelayE=cdDelayE/np.linalg.norm(cdDelayE)

    cdDelayM=np.mean(cdMat[24:32,:],axis=0)
    cdDelayM=cdDelayM/np.linalg.norm(cdDelayM)

    cdDelayL=np.mean(cdMat[32:40,:],axis=0)
    cdDelayL=cdDelayL/np.linalg.norm(cdDelayL)

    (proj1E,proj1M,proj1L)=(S1_mm @ cdDelayE, S1_mm @ cdDelayM,S1_mm @ cdDelayL)
    (proj2E,proj2M,proj2L)=(S2_mm @ cdDelayE, S2_mm @ cdDelayM,S2_mm @ cdDelayL)

    fig=plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(proj1E,proj1M,proj1L,'-r')
    ax.plot(proj2E,proj2M,proj2L,'-b')


    cdDelayE=np.mean(cdMat[16:24,:],axis=0)
    cdDelayE=cdDelayE/np.linalg.norm(cdDelayE)

    cdDelayL=np.mean(cdMat[32:40,:],axis=0)
    cdDelayL=cdDelayL/np.linalg.norm(cdDelayL)

    (proj1E,proj1L)=(S1_mm @ cdDelayE,S1_mm @ cdDelayL)
    (proj2E,proj2L)=(S2_mm @ cdDelayE,S2_mm @ cdDelayL)

    testmark=[40,44] if delay==6 else [28,32]

    fig=plt.figure()
    ax = fig.add_subplot(projection='3d')
    h1=ax.plot(proj1E,proj1L,range(56),'-r')
    h2=ax.plot(proj2E,proj2L,range(56),'-b')

    hs=[ax.plot(proj1E[t],proj1L[t],t,'o',c=[0.5,0.5,0.5]) for t in [12,16]]
    [ax.plot(proj2E[t],proj2L[t],t,'o',c=[0.5,0.5,0.5]) for t in [12,16]]

    ht=[ax.plot(proj1E[t],proj1L[t],t,'ko') for t in testmark]
    [ax.plot(proj2E[t],proj2L[t],t,'ko') for t in testmark]

    ax.set_zticks([12.5, 32.5, 52.5])
    ax.set_zticklabels([0, 5, 10])
    ax.set_zlabel('Time (s)')
    ax.set_xlabel('Early delay CD proj. (a.u.)')
    ax.set_ylabel('Late delay CD proj. (a.u.)')
    plt.legend([h1[0],h2[0],hs[0][0],ht[0][0]], ['S1 trial','S2 trial','Sample','Test'])


def cd_projection(delay=6):
    fr_per_su = util.get_dataset()
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])

    cdMat=(S1_mm-S2_mm)
    cdDelayE=np.mean(cdMat[16:24,:],axis=0)
    cdDelayE=cdDelayE/np.linalg.norm(cdDelayE)

    cdDelayM=np.mean(cdMat[24:32,:],axis=0)
    cdDelayM=cdDelayM/np.linalg.norm(cdDelayM)

    cdDelayL=np.mean(cdMat[32:40,:],axis=0)
    cdDelayL=cdDelayL/np.linalg.norm(cdDelayL)

    (proj1E,proj1M,proj1L)=(S1_mm @ cdDelayE, S1_mm @ cdDelayM,S1_mm @ cdDelayL)
    (proj2E,proj2M,proj2L)=(S2_mm @ cdDelayE, S2_mm @ cdDelayM,S2_mm @ cdDelayL)

    fig=plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(proj1E,proj1M,proj1L,'-r')
    ax.plot(proj2E,proj2M,proj2L,'-b')


    cdDelayE=np.mean(cdMat[16:24,:],axis=0)
    cdDelayE=cdDelayE/np.linalg.norm(cdDelayE)

    cdDelayL=np.mean(cdMat[32:40,:],axis=0)
    cdDelayL=cdDelayL/np.linalg.norm(cdDelayL)

    (proj1E,proj1L)=(S1_mm @ cdDelayE,S1_mm @ cdDelayL)
    (proj2E,proj2L)=(S2_mm @ cdDelayE,S2_mm @ cdDelayL)

    testmark=[40,44] if delay==6 else [28,32]

    fig=plt.figure()
    ax = fig.add_subplot(projection='3d')
    h1=ax.plot(proj1E[8:45],proj1M[8:45],proj1L[8:45],'-r')
    h2=ax.plot(proj2E[8:45],proj2M[8:45],proj2L[8:45],'-b')

    interp=[np.interp(np.array([12,16]+testmark)-0.5,
                      range(oneproj.shape[0]),
                      oneproj)
            for oneproj in [proj1E,proj1M,proj1L,proj2E,proj2M,proj2L]]

    hs=[ax.plot(interp[0][t],interp[1][t],interp[2][t],'o',c=[0.5,0.5,0.5]) for t in [0,1]]
    [ax.plot(interp[3][t],interp[4][t],interp[5][t],'o',c=[0.5,0.5,0.5]) for t in [0,1]]

    ht=[ax.plot(interp[0][t],interp[1][t],interp[2][t],'ko') for t in [2,3]]
    [ax.plot(interp[3][t],interp[4][t],interp[5][t],'ko') for t in [2,3]]

    # ax.set_zticks([12.5, 32.5, 52.5])
    # ax.set_zticklabels([0, 5, 10])
    ax.set_zlabel('Late delay CD proj. (a.u.)')
    ax.set_xlabel('Early delay CD proj. (a.u.)')
    ax.set_ylabel('Mid delay CD proj. (a.u.)')
    plt.legend([h1[0],h2[0],hs[0][0],ht[0][0]], ['S1 trial','S2 trial','Sample','Test'])

def cd_distance(delay=6):
    fr_per_su = util.get_dataset()
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])

    cdMat=(S1_mm-S2_mm)

    cdDelayE=np.mean(cdMat[16:24,:],axis=0)
    cdDelayE=cdDelayE/np.linalg.norm(cdDelayE)

    cdDelayM=np.mean(cdMat[24:32,:],axis=0)
    cdDelayM=cdDelayM/np.linalg.norm(cdDelayM)

    cdDelayL=np.mean(cdMat[32:40,:],axis=0)
    cdDelayL=cdDelayL/np.linalg.norm(cdDelayL)

    (proj1E,proj1M,proj1L)=(S1_mm @ cdDelayE, S1_mm @ cdDelayM,S1_mm @ cdDelayL)
    (proj2E,proj2M,proj2L)=(S2_mm @ cdDelayE, S2_mm @ cdDelayM,S2_mm @ cdDelayL)
    eculidian=[np.linalg.norm([proj1E[t]-proj2E[t],proj1M[t]-proj2M[t],proj1L[t]-proj2L[t]]) for t in range(len(proj1E))]

    fig=plt.figure()
    he=plt.plot(proj1E-proj2E,'-r')
    hm=plt.plot(proj1M-proj2M,'-m')
    hl=plt.plot(proj1L-proj2L,'-b')
    hecu=plt.plot(eculidian,'-k')
    ax=plt.gca()
    ax.set_xticks([12.5, 32.5, 52.5])
    ax.set_xticklabels([0, 5, 10])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('S1-S2 C.D. distance (a.u.)')
    testmark=[40,44] if delay==6 else [28,32]
    [plt.axvline(x,ls='--',c='k') for x in np.array([12,16]+testmark)+0.5]
    plt.legend((he[0],hm[0],hl[0],hecu[0]),('Early CD','Mid CD','Late CD','Euclidian'))