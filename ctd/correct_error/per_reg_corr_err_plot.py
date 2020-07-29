# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt


def gauss_average(x):
    return np.convolve(x, [0.1968, 0.6063, 0.1968], "same")


efstr=np.load('per_region_ctd_6_100_SVC_error.npz')
e_reg=efstr['reg_set'].tolist()
e_mat=efstr['decode_mat']
e_no=efstr['actual_number']


cfstr=np.load('per_region_ctd_6_100_SVC.npz')
c_reg=cfstr['reg_set'].tolist()
c_mat=cfstr['decode_mat']
c_no=cfstr['actual_number']


fullset=list(set(c_reg+e_reg))

plt.close('all')
for i in range(len(fullset)):
    print(i)
    curr_reg=fullset[i]
    
    if (curr_reg in c_reg) and (curr_reg in e_reg):
        (fh,ax)=plt.subplots(1,3,figsize=(6,2),dpi=300)     
        curr_emat=e_mat[e_reg.index(curr_reg),:,:,:]
        curr_cmat=c_mat[c_reg.index(curr_reg),:,:,:]
        
        im=ax[1].imshow(np.mean(curr_emat,axis=0),origin='lower',cmap='jet',vmin=25,vmax=75)
        ax[0].imshow(np.mean(curr_cmat,axis=0),origin='lower',cmap='jet',vmin=25,vmax=75)
        
        for oneax in ax[:2]:
            [oneax.axhline(x, color='k', ls='--') for x in [3.5, 7.5]]
            [oneax.axvline(x, color='k', ls='--') for x in [3.5, 7.5]]
            oneax.set_xticks([3.5,23.5])
            oneax.set_yticks([3.5,23.5])
            oneax.set_xticklabels([0,5])
            oneax.set_yticklabels([0,5])

        ax[1].set_yticks([])
        
        ax[0].set_title('correct')
        ax[1].set_title('error')
        
        mmc=[np.mean(curr_cmat[:,x,x]) for x in range(32)]
        mme=[np.mean(curr_emat[:,x,x]) for x in range(32)]
        
        ax[2].plot(gauss_average(mmc),'r-')
        ax[2].plot(gauss_average(mme),'k-')
        ax[2].set_ylim([40,100])
        [ax[2].axvline(x,color='k',ls='--') for x in [3.5,7.5]]
        ax[2].axhline(50,color='k',ls=':')
        ax[2].set_xticks([3.5,23.5])
        ax[2].set_xticklabels([0,5])
        
        fh.suptitle(curr_reg)
        fh.savefig(f'correct_error_{curr_reg}.pdf',bbox_inches='tight')
        # plt.show()


plt.close('all')
(fh,ax)=plt.subplots(1,1,figsize=(4,4),dpi=300)
for i in range(len(fullset)):
    curr_reg=fullset[i]

    if (curr_reg in c_reg) and (curr_reg in e_reg):

        curr_emat=e_mat[e_reg.index(curr_reg),:,:,:]
        curr_cmat=c_mat[c_reg.index(curr_reg),:,:,:]

        mmc=[np.mean(curr_cmat[:,x,x]) for x in range(8,32)]
        mme=[np.mean(curr_emat[:,x,x]) for x in range(8,32)]

        mmmc=np.mean(mmc)
        mmme=np.mean(mme)

        ax.plot([0,100],[0,100],'k:',lw=0.5)

        ax.plot(mmmc,mmme,'r.')
        ax.text(mmmc,mmme,curr_reg,ha='center')
        ax.set_xlim([50,90])
        ax.set_ylim([40,75])
        ax.set_xlabel('correct trials')
        ax.set_ylabel('error trials')
        ax.set_title('decoding accuracy')
plt.show()



plt.close('all')
(fh,ax)=plt.subplots(1,1,figsize=(4,4),dpi=300)
for i in range(len(fullset)):
    curr_reg=fullset[i]

    if (curr_reg in c_reg) and (curr_reg in e_reg):

        curr_emat=e_mat[e_reg.index(curr_reg),:,8:32,8:32].copy()
        curr_cmat=c_mat[c_reg.index(curr_reg),:,8:32,8:32].copy()

        diagc=[np.mean(curr_cmat[:,x,x]) for x in range(24)]
        diage=[np.mean(curr_emat[:,x,x]) for x in range(24)]

        for dd in range(24):
            curr_cmat[:,dd,dd]=np.nan
            curr_emat[:,dd,dd]=np.nan

        mratioc=np.mean(diagc)/np.nanmean(curr_cmat)
        mratioe=np.mean(diage)/np.nanmean(curr_emat)

        # ax.plot([0,100],[0,100],'k:',lw=0.5)

        ax.plot(mratioc,mratioe,'r.')
        ax.text(mratioc,mratioe,curr_reg,ha='center')
        # ax.set_xlim([50,90])
        # ax.set_ylim([40,75])
        ax.set_xlabel('correct trials')
        ax.set_ylabel('error trials')
        ax.set_title('persistency (diag/offdiag)')
plt.show()
        
        # im=ax[1].imshow(np.mean(curr_emat,axis=0),origin='lower',cmap='jet',vmin=25,vmax=75)
        # ax[0].imshow(np.mean(curr_cmat,axis=0),origin='lower',cmap='jet',vmin=25,vmax=75)
        
        # for oneax in ax[:2]:
        #     [oneax.axhline(x, color='k', ls='--') for x in [3.5, 7.5]]
        #     [oneax.axvline(x, color='k', ls='--') for x in [3.5, 7.5]]
        #     oneax.set_xticks([3.5,23.5])
        #     oneax.set_yticks([3.5,23.5])
        #     oneax.set_xticklabels([0,5])
        #     oneax.set_yticklabels([0,5])

        # ax[1].set_yticks([])
        
        # ax[0].set_title('correct')
        # ax[1].set_title('error')
        
        # mmc=[np.mean(curr_cmat[:,x,x]) for x in range(32)]
        # mme=[np.mean(curr_emat[:,x,x]) for x in range(32)]
        
        # ax[2].plot(gauss_average(mmc),'r-')
        # ax[2].plot(gauss_average(mme),'k-')
        # ax[2].set_ylim([40,100])
        # [ax[2].axvline(x,color='k',ls='--') for x in [3.5,7.5]]
        # ax[2].axhline(50,color='k',ls=':')
        # ax[2].set_xticks([3.5,23.5])
        # ax[2].set_xticklabels([0,5])
        
        # fh.suptitle(curr_reg)
        # fh.savefig(f'correct_error_{curr_reg}.png',bbox_inches='tight')