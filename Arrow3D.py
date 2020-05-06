# -*- coding: utf-8 -*-
"""
Created on Tue May  5 01:10:36 2020

@author: zx
"""


import numpy as np
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

####################################################
# This part is just for reference if
# you are interested where the data is
# coming from
# The plot is at the bottom
#####################################################

# Generate some example data
# mu_vec1 = np.array([0,0,0])
# cov_mat1 = np.array([[1,0,0],[0,1,0],[0,0,1]])
# class1_sample = np.random.multivariate_normal(mu_vec1, cov_mat1, 20)

# mu_vec2 = np.array([1,1,1])
# cov_mat2 = np.array([[1,0,0],[0,1,0],[0,0,1]])
# class2_sample = np.random.multivariate_normal(mu_vec2, cov_mat2, 20)