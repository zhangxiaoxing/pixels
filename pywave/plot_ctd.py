# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 09:48:26 2021

@author: Libra
"""
import matplotlib.pyplot as plt
import numpy as np


def plot_ctd(sus50,trans50,trans1000,delay,repeats):

        (fig, ax) = plt.subplots(1, 3, figsize=[9, 3], dpi=200)
        ax[0].imshow(np.array(sus50).mean(axis=0).mean(axis=0), cmap="jet", aspect="auto", origin='lower', vmin=0,
                     vmax=100)
        ax[0].set_title('50 sustained SU')
        ax[0].set_ylabel('template time (s)')
        ax[1].imshow(np.array(trans50).mean(axis=0).mean(axis=0), cmap="jet", aspect="auto", origin='lower', vmin=0,
                     vmax=100)
        ax[1].set_title('50 transient SU')
        im2 = ax[2].imshow(np.array(trans1000).mean(axis=0).mean(axis=0), cmap="jet", aspect="auto", origin='lower',
                           vmin=0,
                           vmax=100)
        ax[2].set_title('1000 transient SU')
        plt.colorbar(im2, ticks=[0, 50, 100], format="%d")

        for oneax in ax:
            oneax.set_xticks([7.5, 27.5])
            oneax.set_xticklabels([0, 5])
            oneax.set_xlabel('scoring time (s)')
            oneax.set_yticks([7.5, 27.5])
            oneax.set_yticklabels([0, 5])

            if delay == 6:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 35.5, 39.5]]
            elif delay == 3:
                [oneax.axhline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]
                [oneax.axvline(x, color='w', ls=':') for x in [7.5, 11.5, 23.5, 27.5]]

        fig.savefig(f'ctd_{delay}_{repeats}.png', bbox_inches='tight')

        plt.show()