total_n=0;
for sessIdx=1:113
per_samp_trial=1;

% keyboard
if (~exist('path_list','var')) || (exist('to_load','var') && to_load)
    addpath(fullfile('k:','code','jpsth','npy-matlab-master','npy-matlab'))
    addpath('k:\code\jpsth\fieldtrip-20200320')
    ft_defaults
    
    load k:\code\jpsth\114_sorted_file_path.mat
    if isunix
        cd('~/pixels/jpsth')
        homedir='/home/zx/neupix/wyt';
    elseif ispc
        homedir='k:\neupix\wyt';
    end
    sus_trans=h5read('../transient_6.hdf5','/sus_trans'); %export_arr = np.vstack((sust, transient, switched, unclassified, early_in_6s, late_in_6s, prefer_s))
    reg_list=h5read('../transient_6.hdf5','/reg');
    cid_list=h5read('../transient_6.hdf5','/cluster_id');
    path_list=h5read('../transient_6.hdf5','/path');
end



folder=regexp(sorted_fpath{sessIdx},'(\w|\\|-)*','match','once');
stub=regexp(folder,'.*(?=\\)','match','once');
sesssel=startsWith(path_list,stub);
p1sel=contains(path_list,'imec0');
p2sel=contains(path_list,'imec1');
memorysel=any(sus_trans(:,1:2)>0,2);
p1id=cid_list(sesssel & p1sel & memorysel);
reg1=reg_list(sesssel & p1sel & memorysel);
pref1=max(sus_trans(sesssel & p1sel & memorysel,8:13),[],2);
p2id=cid_list(sesssel & p2sel & memorysel)+10000;
reg2=reg_list(sesssel & p2sel & memorysel);
pref2=max(sus_trans(sesssel & p2sel & memorysel,8:13),[],2);

[folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
currmodel='selec';
sustIds=[];nonselIds=[];
transIds=int32([p1id;p2id]);
regs=[reg1;reg2];
prefs=[pref1;pref2];
%% improved pipeline
% greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), regs);
% labeled=~strcmp(deblank(regs),'Unlabeled');
% transIds=transIds(greymatter & labeled);
% regs=regs(greymatter & labeled);
%% compatible pipeline
greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), regs);
transIds=transIds(greymatter);
regs=regs(greymatter);
prefs=prefs(greymatter);

total_n=total_n+nnz(~strcmp(regs,'Unlabeled'));
% msize=numel(transIds);
[avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
delay=6;

TT1 = find(spktrial.trialinfo(:,5)==4 & spktrial.trialinfo(:,8)==delay);
TT2 = find(spktrial.trialinfo(:,5)==8 & spktrial.trialinfo(:,8)==delay);

% TT = [TT1(1:per_samp_trial);TT2(1:per_samp_trial)];
TT=[TT1;TT2];
trialInfo=spktrial.trialinfo(TT,:);
%% reduced number of trials for debug purpose only
%     TT = TT(1);
%%
% TT=32
NN=size(spktrial.label,1);
bin=[1,7];
spiketrains=cell(1,NN);
spiketrials=cell(1,NN);
su_pref=nan(1,NN);
firingrate=nan(NN,length(TT));
for i=1:NN
    sptr=struct();
    sptr.times=[];
    trials=[];
    for t=1:length(TT)
        ttimes=spktrial.time{i}...
            (spktrial.trial{i}==TT(t) ...
            & spktrial.time{i}<bin(2) ...
            & spktrial.time{i}>=bin(1)).*1000+(t-1)*7000;
        ttrials=ones(size(ttimes))*TT(t);
        sptr.times=[sptr.times,ttimes];
        trials=[trials,ttrials];
        firingrate(i,t)=length(ttimes)./diff(bin);
    end
    sptr.times_units='ms';
    sptr.t_start=1000;
    sptr.t_start_units='ms';
    sptr.t_stop=7*length(TT)*1000;
    sptr.t_stop_units='ms';
    spiketrains{i}=sptr;
    spiketrials{i}=trials;
end
seg.spiketrains=spiketrains;
block=struct();
block.name=num2str(sessIdx);
block.segments={seg};
save(sprintf('spktN13_%d.mat',sessIdx),...
    'block','transIds','sessIdx','regs',...
    'spiketrials','trialInfo','firingrate',...
    'prefs','-v7');


end


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
trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')','full');
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
        if strcmp(mode,'full')
            out=facSeq;
        else
            out=facSeq(all(facSeq(:,9:10),2),:);
        end
    else
        out=[];
    end
end
end

