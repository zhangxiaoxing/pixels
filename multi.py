# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 21:56:58 2019

@author: DELL
"""

import pyKilosort

outs=[]
sessions=[]
(out,session)=pyKilosort.runInDir('D:/Data/M35_20191107_g0/M35_20191107_g0_imec0')
outs.append(out)
sessions.append(session)
(out,session)=pyKilosort.runInDir('D:/Data/M35_20191107_g0/M35_20191107_g0_imec0')
outs.append(out)
sessions.append(session)