ed%idces from hebb_pattern_showcase.m
if ~exist('rings','var')
    if isunix
        addpath(fullfile('npy-matlab-master','npy-matlab'))
        addpath('fieldtrip-20200320')
    else
        addpath(fullfile('K:','Lib','npy-matlab-master','npy-matlab'))
        addpath('K:\Lib\fieldtrip-20200320')
    end
    ft_defaults
    load rings.mat
    
    load 114_sorted_file_path.mat
    delay=6;
    if isunix
        cd('~/pixels/jpsth')
        homedir='/home/zx/neupix/wyt';
    elseif ispc
        homedir='k:\neupix\wyt';
    end
end

%%%%%%%%%%%%%%%%%%%%%
% midx=1;
% active=true;
%%%%%%%%%%%%%%%%%%%%
if active
    rings=rings;
    prefix='ring_freq';
else
    rings=rings_inact;
    prefix='ring_freq_inact';
end
%%%%%%%%%%%%%%%%%%%%%

r1=[];
bins1=[];
for bin=1:6
    bindata=cell2mat(rings(midx,:,bin,1)');
    r1=[r1;bindata];
    bins1=[bins1;ones(size(bindata,1),1)*bin];
end
r2=[];
bins2=[];
for bin=1:6
    bindata=cell2mat(rings(midx,:,bin,2)');
    r2=[r2;bindata];
    bins2=[bins2;ones(size(bindata,1),1)*bin];
end
r1(:,midx+3)=1;
r2(:,midx+3)=2;
r=[r1;r2];
bins=[bins1;bins2];
if ~active
    [~,IA,~]=unique(r(:,1:midx+2),'rows'); %merge same ring in different bins
    r=r(IA,:);
    bins=bins(IA,:);
end
ringsm=r(:,1:midx+2);
samps=r(:,midx+3);
all_rings=cell(0,8);
err_rings=cell(0,8);

for ssidx=srange
    if ssidx>size(ringsm,1)
        break;
    end
    sessIdx=idivide(int32(ringsm(ssidx,1)),int32(100000));
    preferedSample=samps(ssidx);
    folder=regexp(sorted_fpath{sessIdx},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
    currmodel='selec';
    sustIds=[];nonselIds=[];
    transIds=int32(rem(ringsm(ssidx,:),100000))';
    
    %% %%%%%% localize %%%%%%%%%%%%%
    freg=load('0116_memory_conn_chain_duo_6s_1_2.mat','pair_reg','pair_chain');
    load reg_keep.mat
    suid=ringsm(ssidx,:);
    reg=[];
    for i1=suid(:)'
        reg=[reg;{i1,reg_set{freg.pair_reg(find(freg.pair_chain==i1,1))}}];
    end
    %% end of localize
    if numel(unique(reg(:,2)))==1
        continue
    end
    
    msize=numel(transIds);
    [avail,spktrial,spktrial_err]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
    if ~avail
        continue
    end
    S1Trials = find(spktrial.trialinfo(:,5)==4 & spktrial.trialinfo(:,8)==delay);
    S2Trials = find(spktrial.trialinfo(:,5)==8 & spktrial.trialinfo(:,8)==delay);
    [S1FreqMat,S1RingSpkCnt]=genFreqMat(S1Trials,msize,spktrial,active,bins(ssidx));
    [S2FreqMat,S2RingSpkCnt]=genFreqMat(S2Trials,msize,spktrial,active,bins(ssidx));
    all_rings(end+1,:)={ssidx,sessIdx,transIds,preferedSample,S1RingSpkCnt,S2RingSpkCnt,S1FreqMat,S2FreqMat};
    
    
    S1Trials = find(spktrial_err.trialinfo(:,5)==4 & spktrial_err.trialinfo(:,8)==delay);
    S2Trials = find(spktrial_err.trialinfo(:,5)==8 & spktrial_err.trialinfo(:,8)==delay);
    [S1FreqMat,S1RingSpkCnt]=genFreqMat(S1Trials,msize,spktrial_err,active,bins(ssidx));
    [S2FreqMat,S2RingSpkCnt]=genFreqMat(S2Trials,msize,spktrial_err,active,bins(ssidx));
    err_rings(end+1,:)={ssidx,sessIdx,transIds,preferedSample,S1RingSpkCnt,S2RingSpkCnt,S1FreqMat,S2FreqMat};
    
end
blame=datetime();
save(sprintf('%s_%d_%d_%d.mat',prefix,midx+2,min(srange),max(srange)),'all_rings','err_rings','blame')
if isunix
    quit(0);
end

function [freqMat,ringSpkCount]=genFreqMat(trials,msize,spktrial,active,bin)
    ringSpkCount=0;
    freqMat=[];
    for tidx=trials(:)'
        ts_id=[];
        for seqid=1:msize
            tsel=spktrial.trial{seqid}==tidx;
            ts=spktrial.time{seqid}(tsel);
            if active
                ts=ts(ts>=bin & ts<bin+1);
            else
                ts=ts(ts>=-2 & ts<7);
            end
            ts(2,:)=seqid;
            ts_id=[ts_id;ts'];
        end
        [~,s]=sort(ts_id(:,1));
        ts_id=ts_id(s,:);
        ts_id_tagged=relax_tag(ts_id,msize);
        ringSpkCount=ringSpkCount+sum(ts_id_tagged(:,3));
        if active
            if ringSpkCount>0 
                onering=nan(1,9);
                onering(bin+3)=sum(ts_id_tagged(:,3));
                freqMat=[freqMat;onering];
            else
                onering=nan(1,9);
                onering(bin+3)=0;
                freqMat=[freqMat;onering];
            end
        end
        if ~active
            if ringSpkCount>0
                freqMat=[freqMat;...
                    histcounts(ts_id_tagged(ts_id_tagged(:,3)==1,1),-2:7)];
            else
                freqMat=[freqMat;zeros(1,9)];
            end
        end
    end
end
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


function [avail,FT_SPIKE,FT_SPIKE_err]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,model)
sps=30000;
trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')',model);
trials_err=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')','error');
if isempty(trials) || isempty(trials_err)
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

cfg.trl=[trials_err(:,1)-3*sps,trials_err(:,1)+11*sps,zeros(size(trials_err,1),1)-3*sps,trials_err];
FT_SPIKE_err=ft_spike_maketrials(cfg,FT_SPIKE);

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

function out=relax_tag(in,msize)
out=in;
out(:,3)=0;
if size(in,1)<msize
    return
end
skiptag=0;
for i=1:size(in,1)
    if i<skiptag
        continue
    end
    curr_su=in(i,2);
    targets=[(curr_su+1):(curr_su+msize-1),curr_su];
    targets(targets>msize)=targets(targets>msize)-msize;
    if isempty(targets)
        continue
    end
    tsseq=[i,in(i,1:2)];
    for t=targets
        rows=tsseq(end,1)+(1:msize*10);
        rows(rows>length(in))=[];
        if isempty(rows)
            break
        end
        didx=find( ...
            in(rows,2)==t ... %post unit
            & in(rows,1)<tsseq(end,2)+0.01 ...
            & in(rows,1)>tsseq(end,2)+0.0005 ... %matching time window, assuming 1kHz
            ,1);
        if isempty(didx)
            break
        else
            tsseq=[tsseq;tsseq(end,1)+didx,in(tsseq(end,1)+didx,1:2)];
        end
    end
    if length(tsseq)<msize+1
        continue
    else
        out(tsseq(:,1),3)=1;
        skiptag=tsseq(2,1);
    end
end
end
