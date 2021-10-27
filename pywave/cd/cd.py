# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:57:04 2021

@author: Libra
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.signal import find_peaks, peak_widths
import scikits.bootstrap as bootstrap

if __package__:
    from pywave import datautils as util
else:
    import sys
    import os
    sys.path.append(os.path.dirname(__file__)+'/..')
    import datautils as util

rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

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
    cdDelayE=np.mean(cdMat[16:20,:],axis=0)
    cdDelayE=cdDelayE/np.linalg.norm(cdDelayE)

    cdDelayM=np.mean(cdMat[26:30,:],axis=0)
    cdDelayM=cdDelayM/np.linalg.norm(cdDelayM)

    cdDelayL=np.mean(cdMat[36:40,:],axis=0)
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
    # fr_error=util.get_dataset(correct_error='error')
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])

    # S1_mm_err=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_error])
    # S2_mm_err=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_error])

    cdMat=(S1_mm-S2_mm)

    cdDelayE=np.mean(cdMat[16:20,:],axis=0)
    cdDelayE=cdDelayE/np.linalg.norm(cdDelayE)

    cdDelayM=np.mean(cdMat[26:30,:],axis=0)
    cdDelayM=cdDelayM/np.linalg.norm(cdDelayM)

    cdDelayL=np.mean(cdMat[36:40,:],axis=0)
    cdDelayL=cdDelayL/np.linalg.norm(cdDelayL)

    (proj1E,proj1M,proj1L)=(S1_mm @ cdDelayE, S1_mm @ cdDelayM,S1_mm @ cdDelayL)
    (proj2E,proj2M,proj2L)=(S2_mm @ cdDelayE, S2_mm @ cdDelayM,S2_mm @ cdDelayL)

    # (proj1Eerr,proj1Merr,proj1Lerr)=(S1_mm_err @ cdDelayE, S1_mm_err @ cdDelayM,S1_mm_err @ cdDelayL)
    # (proj2Eerr,proj2Merr,proj2Lerr)=(S2_mm_err @ cdDelayE, S2_mm_err @ cdDelayM,S2_mm_err @ cdDelayL)

    # eculidian=[np.linalg.norm([proj1E[t]-proj2E[t],proj1M[t]-proj2M[t],proj1L[t]-proj2L[t]]) for t in range(len(proj1E))]

    fig=plt.figure()
    he=plt.plot(proj1E-proj2E,'-r')
    hm=plt.plot(proj1M-proj2M,'-m')
    hl=plt.plot(proj1L-proj2L,'-b')

    # hee=plt.plot(proj1Eerr-proj2Eerr,'--r')
    # hme=plt.plot(proj1Merr-proj2Merr,'--m')
    # hle=plt.plot(proj1Lerr-proj2Lerr,'--b')

    # hecu=plt.plot(eculidian,'-k')
    ax=plt.gca()
    ax.set_xticks([12.5, 32.5, 52.5])
    ax.set_xticklabels([0, 5, 10])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('S1-S2 C.D. distance (a.u.)')
    testmark=[40,44] if delay==6 else [28,32]
    [plt.axvline(x,ls='--',c='k') for x in np.array([12,16]+testmark)+0.5]
    plt.legend((he[0],hm[0],hl[0]),('Early CD','Mid CD','Late CD'))

def cd_FWHM(delay=6):
    fr_per_su = util.get_dataset()
    S1_mm=np.hstack([np.mean(onesu[f'S1_{delay}'],axis=2) for onesu in fr_per_su])
    S2_mm=np.hstack([np.mean(onesu[f'S2_{delay}'],axis=2) for onesu in fr_per_su])
    cdMat=(S1_mm-S2_mm)

    #TODO:traverse bin

    widths=[]
    binws=[8,6,4,3,2,1]
    for binw in binws:
        widthbin=[]
        # fig=plt.figure()
        for onset in range(16,40,binw):
            cdBin=np.mean(cdMat[onset:onset+binw,:],axis=0)
            cdBin=cdBin/np.linalg.norm(cdBin)
            proj1Bin=S1_mm @ cdBin
            proj2Bin=S2_mm @ cdBin
            curve=proj1Bin-proj2Bin
            # plt.plot(curve)
            (peaks,prop)=find_peaks(curve,distance=56)
            (w,_,_,_)=peak_widths(curve, peaks)
            widthbin.append(w)
        widths.append(widthbin)

    ci=np.array([bootstrap.ci(d,np.mean,n_samples=500) for d in widths])*0.25

    plt.figure()
    plt.fill_between(np.array(binws)*0.25,ci[:,0],ci[:,1],color="r", alpha=0.2)
    plt.plot(np.array(binws)*0.25, [np.mean(x)*0.25 for x in widths],'-r')
    plt.plot([0,3],[0,3],':k')
    ax=plt.gca()
    ax.set_xlabel('CD window (s)')
    ax.set_ylabel('FWHM of CD projection (s)')
    ax.set_xlim((0,2))
    ax.set_ylim((0,4))

if __name__=="__main__":
    cd_projection(6);