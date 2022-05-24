# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:48:08 2020

@author: Libra
"""

import matplotlib
import pandas as pd
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.cluster import AgglomerativeClustering


class opgene_tsne:
    def __init__(self):
        self.rois=[]
        self.feat_arr=[]
        self.tsne=[]
        self.colors = ["r", "g", "b", "c", "m"]
        self.assignment=[]


    def calc_tsne(self,feat_arr):
        self.feat_arr=feat_arr
        tsne = TSNE(learning_rate=10, perplexity=5, method="exact", n_iter=10000)
        self.tsne = tsne.fit_transform(self.feat_arr)
        return self.tsne

    def calc_per_su_tsne(self,feat_arr):
        self.feat_arr=feat_arr
        tsne = TSNE()
        self.tsne = tsne.fit_transform(self.feat_arr)
        return self.tsne

    def plot_tsne_agglo(self,n):
        agg = AgglomerativeClustering(n_clusters=n)
        self.assignment = agg.fit_predict(self.tsne)


        fh = plt.figure(figsize=(7.5, 7.5))
        plt.ylim(self.tsne[:, 1].min() - 10, self.tsne[:, 1].max() + 10)
        plt.xlim(self.tsne[:, 0].min() - 10, self.tsne[:, 0].max() + 10)
        for i in range(len(self.rois)):
            plt.plot(self.tsne[i, 0], self.tsne[i, 1], "o", color=self.colors[self.assignment[i]])
            plt.text(self.tsne[i, 0], self.tsne[i, 1] - 5, self.rois[i], ha="center", va="top")
        ax = plt.gca()
        ax.set_title("t-SNE")
        plt.show()
        fh.savefig("tsne_agglo.png", dpi=300, bbox_inches="tight")


    def plot_feature_decomp(self):

        segs = [range(36), range(36, 56), range(60, 88), range(88, 136)]
        clusters = []
        for i in range(5):
            clusters.append(np.mean(self.feat_arr[self.assignment == i, :], axis=0))

        (fh, axes) = plt.subplots(5, 4, figsize=(10, 10), dpi=300)
        for row in range(5):
            for col in range(4):
                axes[row, col].plot(clusters[row][segs[col]], color=self.colors[row])
                axes[row, col].set_ylim((0, 1))
                if col==0 or col==3:
                    yspan = (0, 1)
                    [
                        axes[row,col].plot([x, x], yspan, "--", color="gray")
                        for x in np.array([1, 2, 5, 6]) * 4 - 0.5
                    ]
                    axes[row,col].set_ylim(yspan)
                    axes[row,col].set_xticks([3.5, 7.5, 19.5, 23.5])
                    axes[row,col].set_xticklabels(["S", "+1", "T", "+1"])
                else:
                    yspan = (0, 1)
                    [
                        axes[row,col].plot([x, x], yspan, "--", color="gray")
                        for x in np.array([1, 2]) * 4 - 0.5
                    ]
                    axes[row,col].set_ylim(yspan)
                    axes[row,col].set_xticks([3.5, 7.5])
                    axes[row,col].set_xticklabels(["T", "+1"])

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()
        fh.savefig("match.png", dpi=300, bbox_inches="tight")


    def plot_opgene_mapping(self,perfKey,vgat):
        #%% mappping opto-gene
        if vgat:
            hm=scipy.io.loadmat('K:\\Mapping\\mapping VGAT\\heatMat.mat')['heatMat']
        else:
            hm=scipy.io.loadmat('K:\\Mapping\\mapping effect size NpHR (virus)\effect size\\opgene_ephys_corr.mat')['nphrHeatMat']
        cols=[x.tolist()[0] for x in hm[0,1:]]
        rows=[x.tolist()[0] for x in hm[1:,0]]
        heatMat=pd.DataFrame(data=hm[1:,1:],
                             index=rows,
                             columns=cols)
        heatMat=heatMat.applymap(lambda x:x[0][0])
        cmap = matplotlib.cm.get_cmap('jet')
        breakpoint()
        perf_sel=heatMat.loc[perfKey,:]
        norm = matplotlib.colors.Normalize(-1, 1)

        # heatMat=matFS['heatMat']
        fh = plt.figure(figsize=(7.5, 7.5))
        plt.ylim(self.tsne[:, 1].min() - 10, self.tsne[:, 1].max() + 10)
        plt.xlim(self.tsne[:, 0].min() - 10, self.tsne[:, 0].max() + 10)
        for i in range(len(self.rois)):
            one_roi=self.rois[i]
            color='gray'
            markersize=5
            if one_roi in perf_sel.index:
                color=cmap(norm(perf_sel[one_roi]))
                markersize=10

            plt.plot(self.tsne[i, 0], self.tsne[i, 1], "o", color=color,markersize=markersize)
            plt.text(self.tsne[i, 0], self.tsne[i, 1] - 5, self.rois[i], ha="center", va="top")
        ax = plt.gca()
        ax.set_title("t-SNE, "+perfKey)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        plt.colorbar(sm,ticks=[-1, 0, 1], format="%d")
        plt.show()
        fh.savefig('map_'+perfKey + ".png", dpi=300, bbox_inches="tight")



        #%% map DM-Hit