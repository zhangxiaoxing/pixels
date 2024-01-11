import os, glob, sys, json
import numpy as np


os.chdir(os.path.expanduser("~/npdata_out"))
lrresults = glob.glob("**/time_match_?.json", recursive=True)
all_results=[]
for oneresult in lrresults:
    with open(oneresult, 'r') as file:
        lr=json.load(file)
        all_results.append(lr)
all_results=np.array(all_results)
print('max and mean slope, intercept, r, p, stderr')
print(np.max(all_results,axis=0))
print(np.min(all_results,axis=0))


match_results = glob.glob("**/trials_match.npy", recursive=True)
for one_match in match_results:
    match_arr=np.load(one_match)
    for midx in range(1,match_arr.shape[1]-1):
        reftrial=match_arr.shape[0]
        mismatch=np.sum(np.logical_not(match_arr[:,midx]))
        if mismatch>0 or reftrial<200:
            print((midx,mismatch,reftrial,one_match))
            breakpoint()

