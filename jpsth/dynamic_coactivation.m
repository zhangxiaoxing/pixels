keyboard
sample=1;
% idcesidx=1;
load 114_sorted_file_path.mat
load hebb_pattern.mat
if sample==1
    hebbPattern=hebbPatternS1;
else
    hebbPattern=hebbPatternS2;
end

ssidx=idces(1);
sessIdx=idivide(int32(hebbPattern{ssidx,1}(1)),int32(100000));
hidx=ssidx;

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
%% %%%%%% localize %%%%%%%%%%%%%
freg=load('0831_selec_conn_chain_duo_6s_1_2.mat','pair_reg','pair_chain');
load reg_keep.mat
suid=unique([hebbPattern{hidx,1:3}]);
reg=[];
for i1=suid(:)'
    reg=[reg;{i1,reg_set{freg.pair_reg(find(freg.pair_chain==i1,1))}}];
end
%% end of localize

pairs=rem(unique([hebbPattern{hidx,1}(1:2);...
    hebbPattern{hidx,1}(2:3);...
    hebbPattern{hidx,1}([3 1]);...
    hebbPattern{hidx,2}(1:2);...
    hebbPattern{hidx,2}(2:3);...
    hebbPattern{hidx,2}(3:4);...
    hebbPattern{hidx,2}([4 1]);...
    hebbPattern{hidx,3}(1:2);...
    hebbPattern{hidx,3}(2:3);...
    hebbPattern{hidx,3}(3:4);...
    hebbPattern{hidx,3}([4 1])],'rows'),100000);

hebbpool=unique(pairs(:));

pathstub=regexp(sorted_fpath{sessIdx},'^.*?(?=\\)','match','once');
pathsel=(startsWith(path_list,pathstub) & ismember(cid_list,hebbpool));
auc(pathsel,:);


folder=regexp(sorted_fpath{sessIdx},'(\w|\\|-)*','match','once');
[folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,cell(0),homedir);

sustIds=[];nonselIds=[];
transIds=hebbpool;
currmodel='selec';
[avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix

delay=6;
dt=0.1;
counter=0;
xcs=cell(0,2);
for tt=1:dt:7-dt
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

%%%%%%current case 60bin/trial * 20 trials
coact_ring=zeros(3,1200);
for i=1:1200
    if cross_trial_S1(1,i)>0 && cross_trial_S1(6,i)>0 && cross_trial_S1(7,i)>0
        coact_ring(1,i)=1;
    end
    
    if cross_trial_S1(5,i)>0 && cross_trial_S1(4,i)>0 && cross_trial_S1(8,i)>0 && cross_trial_S1(3,i)>0
        coact_ring(2,i)=2;
    end
    
    if cross_trial_S1(7,i)>0 && cross_trial_S1(2,i)>0 && cross_trial_S1(9,i)>0 && cross_trial_S1(6,i)>0
        coact_ring(3,i)=3;
    end
end




fh=figure('Color','w')
imagesc(cross_trial_S1(:,1:1200),[0,5])
colormap(jet)
hold on
arrayfun(@(x) xline(x,':w'),60:60:1199)
set(gca,'XTick',30:60:1199,'XTickLabel',1:20);
cb=colorbar
cb.Label.String='Coactivation';
ylabel('Pair #');
xlabel('Trial #');
fh.Position(3)=2560;
fh.Position(4)=256;
exportgraphics(fh,'dynamic_coactivation_raw.pdf');
print(fh,'dynamic_coactivation_raw.png','-dpng','-r300');

fhc=figure('Color','w')
imagesc(coact_ring(:,1:1200))
hold on
arrayfun(@(x) xline(x,':w'),60:60:1199)
set(gca,'XTick',30:60:1199,'XTickLabel',1:20,'YTick',1:3);
ylabel('Ring #');
xlabel('Trial #');
fhc.Position(3)=2560;
fhc.Position(4)=256;
cmap=colormap();
cmap(1,:)=[0 0 0];
colormap(cmap);
exportgraphics(fhc,'dynamic_coactivation_ring.pdf');
print(fhc,'dynamic_coactivation_ring.png','-dpng','-r300');


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
