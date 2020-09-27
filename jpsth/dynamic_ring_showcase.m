
keyboard % not a good idea to runthrough
load('rings.mat','rings'); %from ring_list.m
for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
rings3=[];
for bin=1:6
    rbin=cell2mat(rings(1,:,bin,1)');
    rbin(:,4)=bin;
    rbin(:,5)=1;
    rings3=[rings3;rbin];
    rbin=cell2mat(rings(1,:,bin,2)');
    rbin(:,4)=bin;
    rbin(:,5)=2;
    rings3=[rings3;rbin];
end
peak_stats=nan(size(rings3,1),3);
count_stats=nan(size(rings3,1),3);
for i=1:length(rings3)
    bb=rings3(i,4);
    if rings3(i,5)==1
        cc=fstr{bb}.conn_chain_S1;
        pp=fstr{bb}.peaks1;        
        tt=fstr{bb}.totalcount_S1;
    elseif rings3(i,5)==2
        cc=fstr{bb}.conn_chain_S2;
        pp=fstr{bb}.peaks2;
        tt=fstr{bb}.totalcount_S2;
    else
        keyboard
    end
    peak_stats(i,:)=[pp(cc(:,1)==rings3(i,1) & cc(:,2)==rings3(i,2)),...
        pp(cc(:,1)==rings3(i,2) & cc(:,2)==rings3(i,3)),...
        pp(cc(:,1)==rings3(i,3) & cc(:,2)==rings3(i,1))];
    count_stats(i,:)=[tt(cc(:,1)==rings3(i,1) & cc(:,2)==rings3(i,2)),...
        tt(cc(:,1)==rings3(i,2) & cc(:,2)==rings3(i,3)),...
        tt(cc(:,1)==rings3(i,3) & cc(:,2)==rings3(i,1))];
end

[mtt,sidx]=sort(mean(count_stats,2),'descend');
peaks=[mtt,mean(peak_stats(sidx,:),2),sidx,rings3(sidx,4:5)];

keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ssidx=sidx(1);
sessIdx=idivide(int32(rings3(ssidx,1)),int32(100000));

addpath(fullfile('npy-matlab-master','npy-matlab'))
addpath('fieldtrip-20200320')
ft_defaults

if isunix
    cd('~/pixels/jpsth')
    homedir='/home/zx/neupix/wyt';
elseif ispc
    homedir='k:\neupix\wyt';
end
sus_trans=h5read('../transient_6.hdf5','/sus_trans');
path_list=h5read('../transient_6.hdf5','/path');
cid_list=h5read('../transient_6.hdf5','/cluster_id');
auc=h5read('../transient_6.hdf5','/auc');

load 114_sorted_file_path.mat

pathstub=regexp(sorted_fpath{sessIdx},'^.*?(?=\\)','match','once');

ringpool=rem(rings3(ssidx,1:3),100000);
% pathsel=(startsWith(path_list,pathstub) & ismember(cid_list,ringpool));
% auc(pathsel,:);


pairs=[ringpool(1:2);...
    ringpool(2:3);...
    ringpool([3 1])];


folder=regexp(sorted_fpath{sessIdx},'(\w|\\|-)*','match','once');
[folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,cell(0),homedir);

sustIds=[];nonselIds=[];
transIds=ringpool(:);
currmodel='selec';
[avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix

delay=6;
dt=0.1;
counter=0;
xcs=cell(0,2);
for tt=1:dt:2-dt
    counter=counter+1;
    bin_range=[tt,(tt+dt)];
    [xc_s1,xc_s2]=plotxcorr(spktrial,delay,bin_range,'yes');
    xcs(end+1,:)={xc_s1.trial,xc_s2.trial};
end

save('dynamic_coactivation.mat','xcs','xc_s1','xc_s2','pairs')
return

cross_trial_S1=[];
for tid=1:size(xc_s1.trial,1)
    coactMat=nan(size(pairs,1),size(xcs,1));
    for pidx=1:length(pairs)
        preIdx=find(strcmp(xc_s1.label,num2str(pairs(pidx,1))));
        postIdx=find(strcmp(xc_s1.label,num2str(pairs(pidx,2))));
        coactMat(pidx,:)=arrayfun(@(x) sum(xcs{x,1}(tid,preIdx,postIdx,1:4)),1:length(xcs));
    end
    cross_trial_S1=[cross_trial_S1,coactMat];
end

cross_trial_S2=[];
for tid=1:size(xc_s2.trial,1)
    coactMat=nan(size(pairs,1),size(xcs,1));
    for pidx=1:length(pairs)
        preIdx=find(strcmp(xc_s2.label,num2str(pairs(pidx,1))));
        postIdx=find(strcmp(xc_s2.label,num2str(pairs(pidx,2))));
        coactMat(pidx,:)=arrayfun(@(x) sum(xcs{x,1}(tid,preIdx,postIdx,1:4)),1:length(xcs));
    end
    cross_trial_S2=[cross_trial_S2,coactMat];
end


f2=figure()
imagesc(cross_trial_S2)
colorbar()

f1=figure()
imagesc(cross_trial_S1>0)
hold on
arrayfun(@(x) xline(x,'--w'),15.5:15:size(cross_trial_S1,2)-0.5)



figure()
co_ring=sum(cross_trial_S2);
histogram(co_ring)
set(gca,'YScale','log')




return 


function [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list,homedir)
    metaFolder=replace(folder,'\','/');
    metaFolder=fullfile(fullfile(homedir,'DataSum'),metaFolder);
    if isfolder(metaFolder)
        spkFolder=replace(metaFolder,'imec1','imec0');
        file=dir(fullfile(spkFolder,'spike_info.mat'));
        if isempty(file)
            folderType=-1;
            file=[];
            spkFolder=[];
            disp('Error processing file 2-tracks');
            disp(metaFolder);
            error_list(end+1,:)={folderType,metaFolder};
%             pause;
            return
        end
        folderType=2;
    else
        metaFolder=replace(metaFolder,'DataSum','DataSum/singleProbe');
        spkFolder=metaFolder;
        file=dir(fullfile(spkFolder,'spike_times.npy'));
        if isempty(file)
            folderType=-1;
            file=[];
            spkFolder=[];
            disp('Error processing file 1-track');
            disp(metaFolder);
            error_list(end+1,:)={folderType,metaFolder};
            return
        end
        folderType=1;        
    end
end

function [avail,out]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,model)
sps=30000;
trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')',model);
if isempty(trials)
    avail=false;
    out=[];
    return
end

%     trials=double(trials);
%     info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7),trials(:,8)];
if strcmp(model, 'full')
    cluster_ids=[sustIds;transIds;nonselIds];
elseif startsWith(model,'selec')
    cluster_ids=[sustIds;transIds];
elseif startsWith(model,'nonsel')
    cluster_ids=nonselIds;
else
    keyboard
end

%  single-unit candidate

if folderType==1
    spkTS=readNPY(fullfile(spkFolder,'spike_times.npy'));
    spkId=readNPY(fullfile(spkFolder,'spike_clusters.npy'));
elseif folderType==2
    fstr=load(fullfile(spkFolder,'spike_info.mat'));
    spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
    spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);
end
FT_SPIKE=struct();

FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
for i=1:numel(cluster_ids)
    FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
end
%  continuous format F T struct file
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end

function out=clearBadPerf(facSeq, mode)
if strcmp(mode, 'error')
    if length(facSeq)>=40
        errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
        out=facSeq(errorsel,:);
    else
        out=[];
    end
else
    if length(facSeq)>=40
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            good=xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0);
            facSeq(i-39:i,10)=good;
            if nnz(good)>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        out=facSeq(all(facSeq(:,9:10),2),:);
    else
        out=[];
    end
end
end

function [Xc_S1,Xc_S2]=plotxcorr(spikeTrials,delay,bin_range,keeptrials)
% https://www.nature.com/articles/nn799
% A role for inhibition in shaping the temporal flow of information in prefrontal cortex
% Christos Constantinidis, Graham V. Williams & Patricia S. Goldman-Rakic 
% Nature Neuroscience volume 5, pages175-180(2002)
% 
% Neuron, Volume 76
% Functional Microcircuit Recruited during Retrieval of Object Association Memory in Monkey Perirhinal Cortex
% Toshiyuki Hirabayashi, Daigo Takeuchi, Keita Tamura, and Yasushi Miyashita


cfg             = [];
cfg.maxlag      = 0.01; % maximum 100 ms
cfg.binsize     = 0.002; % bins of 2 ms
cfg.outputunit  = 'raw'; % make unit area
cfg.latency     = bin_range; % time bin based on sample onset
cfg.vartriallen = 'no'; % allow variable trial lengths
cfg.debias      = 'no';
cfg.keeptrials  = keeptrials;

cfg.trials      = find(spikeTrials.trialinfo(:,5)==4 & spikeTrials.trialinfo(:,8)==delay);
if numel(cfg.trials)<2
    Xc_S1=[];
else
    cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc_S1 = ft_spike_xcorr(cfg,spikeTrials);
end

cfg.trials      = find(spikeTrials.trialinfo(:,5)==8 & spikeTrials.trialinfo(:,8)==delay);
if numel(cfg.trials)<2
    Xc_S2=[];
else
    cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc_S2 = ft_spike_xcorr(cfg,spikeTrials);
end
% compute the shuffled correlogram
end

