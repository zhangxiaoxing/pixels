# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 22:44:29 2020

@author: Libra
"""

import matplotlib
import pandas as pd
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.cluster import AgglomerativeClustering



from sklearn.manifold import TSNE
tsne=TSNE(learning_rate=10, perplexity=7, n_iter=100000);
for i in range(1):
    print('=============',i,'==================')
    tsne_assign=tsne.fit_transform(feat_per_su_arr[i:23697:5,:])
    print(tsne.n_iter_)
    np.save('tsne_assign'+('%02d' % i)+'.npy',tsne_assign)


tsne_list=[]

for i in range(1):
    print('=============',i,'==================')
    tsne_list.append(np.load('tsne_assign'+('%02d' % i)+'.npy'))
    plt.figure(figsize=(5,5),dpi=300)
    plt.plot(tsne_list[i][:,0],tsne_list[i][:,1],'k.',markersize=1,alpha=0.2)
    plt.show()


n_cluster=40
agg = AgglomerativeClustering(n_clusters=n_cluster)
assignment = agg.fit_predict(tsne_list[0])



fh = plt.figure(figsize=(7.5, 7.5))
norm=matplotlib.colors.Normalize(vmin=0,vmax=n_cluster-1)
cmap = matplotlib.cm.get_cmap('jet')
for i in range(n_cluster):
    plt.plot(tsne_assign[assignment==i,0],tsne_assign[assignment==i,1],'.',color=cmap(norm(i)))
ax = fh.gca()
ax.set_title("t-SNE")
plt.show()
fh.savefig("tsne_agglo_per_neuron_1_in_5.png", dpi=300, bbox_inches="tight")




