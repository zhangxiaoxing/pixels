# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 16:49:02 2020

@author: Libra
"""

""" 
    This tutorial shows how to create and render a brainrender scene with some brain regions
"""
import brainrender
brainrender.SHADER_STYLE = 'cartoon'
brainrender.WHOLE_SCREEN = True
from brainrender.scene import Scene
from matplotlib import cm
import h5py
import numpy as np
import os
from brainrender.animation.video import BasicVideoMaker



with h5py.File(os.path.join(r"K:\code", "transient_6.hdf5"), "r") as ffr:
    sus_trans = np.array(ffr["sus_trans"], dtype="int8")
    reg = [x.decode() for x in ffr["reg"]]

reg_set=list(set(reg))
sums=[]
for one_reg in reg_set:
    reg_sel=[x==one_reg for x in reg]
    count=np.sum(reg_sel)
    trans=np.sum(np.logical_and(reg_sel,sus_trans[1,:]))
    sums.append([one_reg,count,trans,trans/count])

dist=[x[3] for x in sums if x[1]>=100]
dmin=np.floor(np.min(dist)*100).astype(np.int)
dmax=np.floor(np.max(dist)*100).astype(np.int)

screenshot_params = dict(
    folder = './screenshots',
    name='br0',
)

# Create a scene
scene = Scene(screenshot_kwargs=screenshot_params,camera="sagittal")
jmap=cm.get_cmap('jet',dmax-dmin+1)
jetmap=jmap(range(dmax-dmin+1))
# Add the whole thalamus in gray
for s in sums:
    if s[1]>=100 and s[0]!='Unlabeled' and s[0]!='root':
        cmIdx=np.floor(s[3]*100).astype(np.int)-dmin
        scene.add_brain_regions([s[0]], colors=jetmap[cmIdx][:3],use_original_color=False, alpha=0.25,add_labels=True)

vmkr = BasicVideoMaker(scene)
vmkr.make_video(azimuth=2, niters=90, duration=9, save_name="testAZ")

print('done')