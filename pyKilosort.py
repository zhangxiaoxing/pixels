# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:45:30 2019

@author: Libra
"""

import shlex
import subprocess
import os

def run(command):
    try:
        result=subprocess.check_output(shlex.split(command),stderr=subprocess.STDOUT)
        return 0,result
    except subprocess.CalledProcessError as e:
        return e.returncode, e.output

def runInDir(path):
    os.chdir(path)
    status, out=run('matlab -noFigureWindows -batch "lwd=pwd();run D:\code\zxSort.m"')
    if status==0:
        cwd=os.getcwd()
        cleanDir=cwd+'_cleaned'
        os.chdir(cleanDir)
        import sys
        sys.path.insert(1,'D:/code/')
        import sync
        import zxPhy
        import parseDPAFR
        
        sync.runsync()
        zxPhy.runPhy()
        parseDPAFR.runParse()    



    
if __name__=="__main__":
    status, out=run('matlab -noFigureWindows -batch "lwd=pwd();run D:\code\zxSort.m"')
    if status==0:
        cwd=os.getcwd()
        cleanDir=cwd+'_cleaned'
        os.chdir(cleanDir)
        import sys
        sys.path.insert(1,'D:/code/')
        import sync
        import zxPhy
        import parseDPAFR
        
        sync.runsync()
        zxPhy.runPhy()
        parseDPAFR.runParse()