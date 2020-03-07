function CQ_transienti(ubegin,uend)
addpath('CQDecoding')
datapath='su_trials_fr_6.hdf5';
count=h5read(datapath,'/count');
transient=zeros(1,count);
for i=ubegin:uend
    tic
    disp(i);
    fileIdx=num2str(i-1);
    S1_trials=h5read(datapath,'/0/S1');
    S2_trials=h5read(datapath,'/0/S2');
    [~,~,transient(i),~,~,~]=TestSecondSelectivityChange_lite(S1_trials,S2_trials,1,2000,1:6,1:6);
    save(['transient_6_',num2str(ubegin),'_',num2str(uend),'.mat'],'transient','i')
    toc
end
end
