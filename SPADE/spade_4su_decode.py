# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 10:34:34 2020

@author: Libra
"""
import numpy as np
import sys
import scipy.stats as stats
import scipy.io
import os.path
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pickle
import scikits.bootstrap as boot
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import MinMaxScaler,StandardScaler,RobustScaler
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def gen_data(trial_thres=50, error_thres=0):
    congru_stp=[]
    incongru_stp=[]
    
    for sess_id in range(114):
        if not os.path.isfile(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id)):
            continue
    
        filt_fpath=r'K:\code\SPADE\results\{}\winlen4\filtered_patterns.npy'.format(sess_id)
        if not os.path.isfile(filt_fpath):
            continue
    
        mat=scipy.io.loadmat(r'K:\code\SPADE\spkt\spktO17_{}.mat'.format(sess_id))
        pref=np.array(mat['prefs'])[:,0]
        rego=mat['regs']
        regs=[x[0][0] for x in rego]
        
        reg_set=np.unique(regs)
        reg_set=reg_set[reg_set!='Unlabeled']
    
        trialInfo=mat['trialInfo']
        
        s1sel=trialInfo[:,4]==4
        s2sel=trialInfo[:,4]==8
        
        goodsel=np.logical_and(trialInfo[:,8],trialInfo[:,9])
        errorsel=np.logical_not(trialInfo[:,9])
        
        s1c_t=np.nonzero(np.logical_and(s1sel,goodsel))[0]
        s2c_t=np.nonzero(np.logical_and(s2sel,goodsel))[0]
        s1e_t=np.nonzero(np.logical_and(s1sel,errorsel))[0]
        s2e_t=np.nonzero(np.logical_and(s2sel,errorsel))[0]
            
        
        
        
        if s1c_t.shape[0]<trial_thres or s2c_t.shape[0]<trial_thres \
            or s1e_t.shape[0]<error_thres or s2e_t.shape[0]<error_thres:
            continue
        
        print([s1c_t.shape[0],s2c_t.shape[0],s1e_t.shape[0],s2e_t.shape[0],])
        
        r=np.load(filt_fpath,allow_pickle=True)
        patterns=r[0]
        
        for idx, patt in enumerate(patterns):
            reg_grp=[regs[x] for x in patt['neurons']]
            #select only inter-region connection
            if np.unique(reg_grp).shape[0]<2:
                # print('local')
                continue
            if np.isin('Unlabeled',reg_grp):
                # print('Unlabeled')
                continue
            
            pref_grp=[pref[x] for x in patt['neurons']]
            
    ### Actual data accessed after this line
    
            times=patt['times']
            hist,edge=np.histogram(times,np.arange(0,7000*(len(trialInfo))+1,1000))
            s1corr=[]
            s1err=[]
            s2corr=[]
            s2err=[]
            
            for t in s1c_t:
                tt=t*7;
                s1corr.append(hist[tt:tt+6])
            for t in s2c_t:
                tt=t*7;
                s2corr.append(hist[tt:tt+6])
            for t in s1e_t:
                tt=t*7;
                s1err.append(hist[tt:tt+6])
            for t in s2e_t:
                tt=t*7;
                s2err.append(hist[tt:tt+6])
            
            if np.unique(pref_grp).shape[0]<2:
                congru_stp.append({'s1corr':s1corr,'s2corr':s2corr,
                                   's1err':s1err,'s2err':s2err})
                # prefered_samp=np.unique(pref_grp)[0]
            else:
                incongru_stp.append({'s1corr':s1corr,'s2corr':s2corr,
                       's1err':s1err,'s2err':s2err})
                # prefered_samp=0
            
    pickle.dump({'congru_stp':congru_stp,'incongru_stp':incongru_stp},
                open(f'stp_decoding_{trial_thres}_{error_thres}.p','wb'))

def load_data(trial_thres=50, error_thres=0):
    return pickle.load(open(f'stp_decoding_{trial_thres}_{error_thres}.p','rb'))
    

# def plot():




def gen_mat(data,STP_n=100,corr_trial=50,err_trial=0):
    if STP_n>0:
        nn=np.random.choice(data,STP_n,replace=False)
    else:
        nn=data
     
    out=[{
          's1corr':[n['s1corr'][x] for x in np.random.choice(len(n['s1corr']),corr_trial,replace=False)],
        's2corr':[n['s2corr'][x] for x in np.random.choice(len(n['s2corr']),corr_trial,replace=False)],
        's1err':[n['s1err'][x] for x in np.random.choice(len(n['s1err']),err_trial,replace=False)],
        's2err':[n['s2err'][x] for x in np.random.choice(len(n['s2err']),err_trial,replace=False)],
      } for n in nn]
     
    return out


def decode(dec_data,corr_err=False):
    rng = np.random.default_rng()
    dec=[]
    shuf_dec=[]
    err_dec=[]
    shuf_err_dec=[]
    scaler = MinMaxScaler()
    trials=len(dec_data[0]['s1corr'])
    kf=KFold(n_splits=trials)
    # clf = SVC(C=0.001,kernel='linear')
    clf=LinearDiscriminantAnalysis()
    for t in [0,2,4]:#range(6):
        X1=np.sum(np.array([n['s1corr'] for n in dec_data])[:,:,t:t+2],axis=2).transpose()
        X2=np.sum(np.array([n['s2corr'] for n in dec_data])[:,:,t:t+2],axis=2).transpose()
        X1err=np.sum(np.array([n['s1err'] for n in dec_data])[:,:,t:t+2],axis=2).transpose()
        X2err=np.sum(np.array([n['s2err'] for n in dec_data])[:,:,t:t+2],axis=2).transpose()
        
        scaler = scaler.fit(np.vstack((X1, X2, X1err, X2err)))


        XE=scaler.transform(np.vstack((X1err, X2err)))
        yE=np.hstack((np.zeros(X1err.shape[0]),np.ones(X2err.shape[0])))
        deckf=[]
        errkf=[]
        shufkf=[]
        for train_idx, test_idx in  kf.split(X1):
            X=scaler.transform(np.vstack((X1[train_idx], X2[train_idx])))
            y=np.hstack((np.zeros(X1.shape[0])[train_idx],np.ones(X2.shape[0])[train_idx]))
            clf.fit(X,y)            
            XT=scaler.transform(np.vstack((X1[test_idx], X2[test_idx])))
            yT=np.hstack((np.zeros(X1.shape[0])[test_idx],np.ones(X2.shape[0])[test_idx]))
            deckf.extend(clf.predict(XT)==yT)
            
            errkf.extend(clf.predict(XE)==yE)
            yS=y.copy()
            rng.shuffle(yS)
            clf.fit(X,yS)
            shufkf.extend(clf.predict(XT)==yT)
        
        dec.append(deckf)
        err_dec.append(errkf)
        shuf_dec.append(shufkf)

            
    return (dec,shuf_dec,err_dec,shuf_err_dec)
        
def plot(data_arr,data_shuf_arr,data_err_arr,file_desc='sample'):
    
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    
    cmm=np.mean(data_arr,axis=1)
    cshufmm=np.mean(data_shuf_arr,axis=1)
    cerrmm=np.mean(data_err_arr,axis=1)
    
    data_boot=[boot.ci(data_arr[b,:], np.mean,n_samples=1000) for b in range(3)]
    data_shuf_boot=[boot.ci(data_shuf_arr[b,:], np.mean,n_samples=1000) for b in range(3)]
    data_err_boot=[boot.ci(data_err_arr[b,:], np.mean,n_samples=1000) for b in range(3)]
    
    (fh, ax) = plt.subplots(1, 1, figsize=(5 / 2.54, 5 / 2.54), dpi=300)
    ax.fill_between(np.arange(1,7,2),[x[0] for x in data_boot],[x[1] for x in data_boot],color='r',alpha=0.2)
    ax.fill_between(np.arange(1,7,2),[x[0] for x in data_shuf_boot],[x[1] for x in data_shuf_boot],color='k',alpha=0.2)
    ax.fill_between(np.arange(1,7,2),[x[0] for x in data_err_boot],[x[1] for x in data_err_boot],color='b',alpha=0.2)
    
    ax.plot(np.arange(1,7,2),np.mean(data_arr,axis=1),'-r')
    ax.plot(np.arange(1,7,2),np.mean(data_shuf_arr,axis=1),'-k')
    ax.plot(np.arange(1,7,2),np.mean(data_err_arr,axis=1),'-b')
    ax.set_xlim((0,6.5))
    ax.set_ylim((0.2,0.8))
    ax.set_xticks((0,5))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Classification accuracy')
    fh.savefig('4su_stp_decoding_{}.pdf'.format(file_desc),bbox_inches='tight')
    return fh

def gen_dec_arr(STP_n=100,rpt=100,trial_thres=50, error_thres=0,corr_err=False):
    data=load_data(trial_thres,error_thres)
    congru_dec=[]
    congru_shuf=[]
    congru_err=[]
    congru_shuf_err=[]
    
    # incongru_dec=[]
    # incongru_shuf=[]
    # incongru_err=[]
    # incongru_shuf_err=[]
    
    for i in range(rpt):
        print(f'repeat {i}')
        congru_data=gen_mat(data['congru_stp'],STP_n=STP_n,corr_trial=trial_thres,err_trial=error_thres)
        (dec,shuf_dec,err_dec,shuf_err_dec)=decode(congru_data,corr_err)
        congru_dec.append(dec)
        congru_shuf.append(shuf_dec)
        congru_err.append(err_dec)
        congru_shuf_err.append(shuf_err_dec)
        
        # incongru_data=gen_mat(data['incongru_stp'],STP_n=STP_n)
        # (inc_dec,inc_shuf_dec,inc_err_dec,inc_shuf_err_dec)=decode(incongru_data)
        # incongru_dec.append(inc_dec)
        # incongru_shuf.append(inc_shuf_dec)
        # incongru_err.append(inc_err_dec)
        # incongru_shuf_err.append(inc_shuf_err_dec)
    
    cong_arr=np.hstack(tuple(congru_dec))
    cong_shuf_arr=np.hstack(tuple(congru_shuf))
    cong_err_arr=np.hstack(tuple(congru_err))
    cong_shuf_err_arr=np.hstack(tuple(congru_shuf_err))
    # incong_arr=np.hstack(tuple(incongru_dec))
    # incong_shuf_arr=np.hstack(tuple(incongru_shuf))
    # incong_err_arr=np.hstack(tuple(incongru_err))
    # incong_shuf_err_arr=np.hstack(tuple(incongru_shuf_err))
    
    incong_arr=None
    incong_shuf_arr=None
    incong_err_arr=None
    incong_shuf_err_arr=None
    
    pickle.dump({'congru':cong_arr,
                 'congru_shuf':cong_shuf_arr,
                 'congru_err':cong_err_arr,
                 'congru_shuf_err':cong_shuf_err_arr,
                 'incongru':incong_arr,
                 'incongru_shuf':incong_shuf_arr,
                 'incongru_err':incong_err_arr,
                 'incongru_shuf_err':incong_shuf_err_arr,},
                open('stp_decoding_{}spt_{}rpt_{}_{}.p'
                     .format(STP_n,rpt,trial_thres,error_thres),'wb'))


### ########
if __name__ == "__main__":
    if True: # Sample decoding
        trial_thres=45
        error_thres=5
        rpt=25
        STP_n=0
        if False:
            gen_data(trial_thres=trial_thres,error_thres=error_thres)
        
        gen_dec_arr(STP_n,rpt=rpt,trial_thres=trial_thres,error_thres=error_thres)
        
        data=pickle.load(open('stp_decoding_{}spt_{}rpt_{}_{}.p'
                              .format(STP_n,rpt,trial_thres,error_thres),'rb'))
        fh=plot(data['congru'],data['congru_shuf'],data['congru_err'],file_desc='congru_sample'.format(STP_n))    
    # if False: #correct/error decoding
    #     trial_thres=12
    #     error_thres=12
    #     rpt=2
    #     STP_n=0
    #     if True:
    #         gen_data(trial_thres=trial_thres,error_thres=error_thres)
    #         # sys.exit()
    #     gen_dec_arr(STP_n,
    #                 rpt=rpt,
    #                 trial_thres=trial_thres,
    #                 error_thres=error_thres,
    #                 corr_err=True)
        
    #     data=pickle.load(open('stp_decoding_{}spt_{}rpt_{}_{}.p'
    #                           .format(STP_n,rpt,trial_thres,error_thres),'rb'))
    #     fh=plot(data['congru_err'],data['congru_shuf_err'],None,file_desc='congru_err_{}'.format(STP_n))