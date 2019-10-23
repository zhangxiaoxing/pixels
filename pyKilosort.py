# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:45:30 2019

@author: Libra
"""

import shlex
import subprocess

def run(command):
    try:
        result=subprocess.check_output(shlex.split(command),stderr=subprocess.STDOUT)
        return 0,result
    except subprocess.CalledProcessError as e:
        return e.returncode, e.output
    
    
if __name__=="__main__":
    status, out=run('matlab -batch "run K:\code\zxSort.m"')
if status!=0:
    import sys
    sys.path.insert(1,'K:/code/')
    import sync
    import zxPhy
    import parseDPAFR
    
    sync.runsync()
    zxPhy.runPhy()
    parseDPAFR.runParse()