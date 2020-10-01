%  A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
%  chan2. Thus, a negative lag represents a spike in the second dimension of
%  stat.xcorr before the channel in the third dimension of stat.stat.
addpath(fullfile('npy-matlab-master','npy-matlab'))
addpath('fieldtrip-20200320')
ft_defaults


currmodel='full';%selec, full, nonsel
prefix='0831';
delay=6;
if ~exist('fsum','var')
    fsum=load(sprintf('0831_selec_XCORR_duo_sums_delay_6_%d_%d_2msbin.mat',bin_range(1),bin_range(2)));
    [~,fidx]=sort(fsum.sums(:,2));
    fsum.sums=fsum.sums(fidx,:);
end

if isunix
    cd('~/pixels/jpsth')
    homedir='/home/zx/neupix/wyt';
    addpath('~/refcode/buzcode/io')
    addpath('~/refcode/buzcode/utilities/')
    addpath('~/refcode/buzcode/analysis/spikes/correlation/')
    addpath('~/refcode/buzcode/analysis/spikes/functionalConnectionIdentification/')
    addpath('~/refcode/buzcode/visualization/')
    addpath('~/refcode/buzcode/externalPackages/FMAToolbox/General/');
    addpath('~/refcode/buzcode/externalPackages/FMAToolbox/Helpers/');

elseif ispc
    homedir='k:\neupix\wyt';
    addpath('K:\refcode\buzcode\io')
    addpath('K:\refcode\buzcode\utilities\')
    addpath('K:\refcode\buzcode\analysis\spikes\correlation\')
    addpath('K:\refcode\buzcode\visualization\')
end

addpath(fullfile('npy-matlab-master','npy-matlab'))
error_list=cell(0);


for i=12:length(fsum.sums)
    if exist('xcorrpause','file')
        disp('paused by file')
        disp(i)
        keyboard
    end

    folder=fsum.sums{i,2};
    [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list,homedir);


    sustIds=fsum.sums{i,3};
    transIds=fsum.sums{i,4};
    nonselIds=[];

    
    if numel(fsum.sums{i,5}.cfg.trials)<20 ||  numel(fsum.sums{i,7}.cfg.trials)<20
        skipflag=true;
        sumsbz={i,folder,sustIds,transIds,[],[],nonselIds}; %per folder save
    else
        skipflag=false;
        [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
        [mono_rez_S1,mono_rez_S2]=sortSpikeIDz(spktrial,delay,bin_range);
        sumsbz={i,folder,sustIds,transIds,mono_rez_S1,mono_rez_S2,nonselIds}; %per folder save
    end
    [overlapS1,overlapS2]=overlap_calc(mono_rez_S1,mono_rez_S2,i,bin_range);
    save(sprintf('%s_%s_BZ_XCORR_duo_f%d_delay_%d_%d_%d_2msbin.mat',prefix,currmodel,i,delay,bin_range(1),bin_range(2)),'sumsbz','-v7.3','overlapS1','overlapS2','folder') %prefix
 	fprintf('%d of %d\n',i,length(fsum.sums))
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

function [mono_rez_S1,mono_rez_S2]=sortSpikeIDz(spktrial,delay,bin_range)
    %% s1
    spikeIDz_S1=[];
    TS_S1=[];
    spikeIDz_S2=[];
    TS_S2=[];
    s1trl = find(spktrial.trialinfo(:,5)==4 & spktrial.trialinfo(:,8)==delay);
    s2trl = find(spktrial.trialinfo(:,5)==8 & spktrial.trialinfo(:,8)==delay);
    for i=1:numel(spktrial.label)
        binsel=spktrial.time{i}>=min(bin_range) & spktrial.time{i}<max(bin_range);
        lbl=str2double(spktrial.label{i});
        trialsel=ismember(spktrial.trial{i},s1trl);
        spkts=(double(spktrial.timestamp{i}(trialsel & binsel))./30000)';
        TS_S1=[TS_S1;spkts];
        spikeIDz_S1=[spikeIDz_S1;repmat([0,lbl,lbl],length(spkts),1)];
        
        trialsel=ismember(spktrial.trial{i},s2trl);
        spkts=(double(spktrial.timestamp{i}(trialsel & binsel))./30000)';
        TS_S2=[TS_S2;spkts];
        spikeIDz_S2=[spikeIDz_S2;repmat([0,lbl,lbl],length(spkts),1)];
    end
    [LUT,UB,UC]=unique(spikeIDz_S1(:,3));
    spikeIDz_S1(:,3)=UC;
    mono_rez_S1 = zx_MonoSynConvClick (spikeIDz_S1,TS_S1,'alpha',0.05);
    mono_rez_S1.LUT=LUT;
    
    
    [LUT,UB,UC]=unique(spikeIDz_S2(:,3));
    spikeIDz_S2(:,3)=UC;
    mono_rez_S2 = zx_MonoSynConvClick (spikeIDz_S2,TS_S2,'alpha',0.05);
    mono_rez_S2.LUT=LUT;
    
end


function [overlapS1,overlapS2]=overlap_calc(mono_rez_S1,mono_rez_S2,sess,bin_range)
    fcon=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin_range(1),bin_range(2)));
    overlapS1=[0,0,size(mono_rez_S1.sig_con,1)];
    for i=1:size(mono_rez_S1.sig_con,1)
        t=mono_rez_S1.LUT(mono_rez_S1.sig_con(i,:))+sess*100000;
%         disp(t);
        overlapS1(1)=overlapS1(1)+any(fcon.conn_chain_S1(:,1)==t(1) & fcon.conn_chain_S1(:,2)==t(2));
        overlapS1(2)=overlapS1(2)+any(fcon.pair_chain(:,1)==t(1) & fcon.pair_chain(:,2)==t(2));
    end
    
    overlapS2=[0,0,size(mono_rez_S2.sig_con,1)];
    for i=1:size(mono_rez_S2.sig_con,1)
        t=mono_rez_S2.LUT(mono_rez_S2.sig_con(i,:))+sess*100000;
%         disp(t);
        overlapS2(1)=overlapS2(1)+any(fcon.conn_chain_S2(:,1)==t(1) & fcon.conn_chain_S2(:,2)==t(2));
        overlapS2(2)=overlapS2(2)+any(fcon.pair_chain(:,1)==t(1) & fcon.pair_chain(:,2)==t(2));
    end
end
