function CQ_transient(nPool)
addpath('CQDecoding')
datapath='su_trials_fr_6.hdf5';
count=h5read(datapath,'/count');
%parpool(nPool);

batchLen=ceil(count/nPool);
batchStart=1:batchLen:count;

futures=parallel.FevalFuture.empty(0,nPool);
for batchIdx=1:nPool
    batchEnd=batchStart(batchIdx)+batchLen-1;
    if batchEnd>count
        batchEnd=count;
    end
    futures(batchIdx) = parfeval(@onebatch,0,batchStart(batchIdx),batchEnd);
end
for batchIdx=1:nPool
    completedIdx = fetchNext(futures);
    fprintf('Got result with index: %d.\n', completedIdx);
end

end



function onebatch(ubegin,uend)
addpath('CQDecoding')
datapath='su_trials_fr_6.hdf5';
count=h5read(datapath,'/count');
transient=zeros(1,count);
for i=ubegin:uend
    tic
    disp(i);
    S1_trials=h5read(datapath,['/',num2str(i-1),'/S1']);
    S2_trials=h5read(datapath,['/',num2str(i-1),'/S2']);
    [~,~,transient(i),~,~,~]=TestSecondSelectivityChange_lite(S1_trials,S2_trials,2000,1:3,1:3);
    save(['transient_6_early3',num2str(ubegin),'_',num2str(uend),'.mat'],'transient','i','ubegin','uend')
    toc
end
end




function [ZStatistics,ShuffledZStatistics,IsTransientDelayFR1,P,IsSignificant1,IsSignificantLess1]=TestSecondSelectivityChange_lite(Samp1TrialsFR...
    ,Samp2TrialsFR,ShuffleTimes,Permuted_BinID,DelayBinID)

Samp1TrialsFR = Samp1TrialsFR(:,Permuted_BinID);
Samp2TrialsFR = Samp2TrialsFR(:,Permuted_BinID);
%% construct difference score matrix for each trial
% Stable and dynamic coding for working memory in primate prefrontal cortex,
% Eelke Spaak,Kei Watanabe,Shintaro Funahashi,and Mark G. Stokes
% construct the FR difference score matrix for each trial and each bin during delay period: nTrialNum X DelayBinNum X DelayBinNum
BinNum=size(Samp1TrialsFR,2);
Sam1DiffScore=zeros(size(Samp1TrialsFR,1),BinNum,BinNum);
for i=1:size(Samp1TrialsFR,1)%go through each trial
    tempTrialBinnedFR=Samp1TrialsFR(i,:);
    for iBin=1:BinNum
        Sam1DiffScore(i,iBin,:)=tempTrialBinnedFR-tempTrialBinnedFR(iBin);
    end
end
Sam2DiffScore=zeros(size(Samp2TrialsFR,1),BinNum,BinNum);
for i=1:size(Samp2TrialsFR,1)%go through each trial
    tempTrialBinnedFR=Samp2TrialsFR(i,:);
    for iBin=1:BinNum
        Sam2DiffScore(i,iBin,:)=tempTrialBinnedFR-tempTrialBinnedFR(iBin);
    end
end
%% construct the Z-stastics (based on ranksum test) for each nTimeBins X nTimeBins matrix
ZStatistics=zeros(BinNum,BinNum);
for i=1:BinNum
    for j=1:BinNum
        if j~=i
            [~,~,stats] =ranksum(Sam1DiffScore(:,i,j),Sam2DiffScore(:,i,j),'method','approximate');
            if isfield(stats,'zval') && ~isnan(stats.zval)
                ZStatistics(i,j)=stats.zval;
            else
                ZStatistics(i,j)=0; 
            end
        end
    end
end
%% construct the shuffled difference score and Z-statistics by permuting the time labels for ShuffleTimes times

ShuffledZStatistics=zeros(ShuffleTimes,BinNum,BinNum);
for iShuffleTimes=1:ShuffleTimes
    tempShuffledZStatistics=ComputeFRDiffScore(Samp1TrialsFR,Samp2TrialsFR,BinNum);
    ShuffledZStatistics(iShuffleTimes,:,:)= tempShuffledZStatistics;
end


DelayBinIndex=find(Permuted_BinID>=min(DelayBinID)&Permuted_BinID<=max(DelayBinID));
%% perform permutation test
[IsSignificant1,IsSignificantLess1,P]=PermutationTest(ZStatistics,ShuffledZStatistics,1);
DelayIsSignificant=IsSignificant1(DelayBinIndex,DelayBinIndex);
DelayIsSignificantLess=IsSignificantLess1(DelayBinIndex,DelayBinIndex);
if sum(sum(DelayIsSignificant+DelayIsSignificantLess))-sum(diag(DelayIsSignificant))-sum(diag(DelayIsSignificantLess))>0
    IsTransientDelayFR1=1;
else
    IsTransientDelayFR1=0;
end
end
