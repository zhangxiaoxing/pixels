addpath('npy-matlab-master\npy-matlab')
addpath('fieldtrip-20200320')
ft_defaults
tfs=importdata('transient_6.csv');
sufs=importdata('su_list.csv',',');

sust=find(tfs.data(1,:)>0);
trans=find(tfs.data(2,:)>0);
counter=[];
done=[];
sums=cell(0,5);
for i=1:length(sust)
    if ismember(sust(i),done)
        continue
    end
    folder=sufs.textdata{sust(i)};
    wffile=fullfile(replace(folder,'D:','D:'),'wf_stats.hdf5');
    if isfile(wffile)
        sustIds=sufs.data(strcmp(sufs.textdata,folder)' & tfs.data(1,:)>0);
        sameFolder=find(strcmp(sufs.textdata,folder)' & tfs.data(1,:)>0);
        done=[done,sameFolder];
        sustCount=numel(sustIds);
        transIds=sufs.data(strcmp(sufs.textdata,folder)' & tfs.data(2,:)>0);
        transCount=numel(transIds);
        if transCount<1
            continue
        end
        [avail,spktrial]=pre_process(replace(folder,'D:','D:'),sustIds,transIds);
        if avail
            xc=plotxcorr(spktrial);
        end
        % waveform data
        
        wfstats=h5read(wffile,'/wf');
        for lblidx=1:size(xc.label,1)
            wfidx=find(wfstats(:,1)==str2double(xc.label{lblidx,1}));
            if ~isempty(wfidx)
                xc.label{lblidx,2}=wfstats(wfidx,:);
                raw_wf_file=fullfile(replace(folder,'D:\neupix\DataSum','D:\neupix\WF\neuropixel'),'waveform.mat');
                raw_fstr=load(raw_wf_file);
                raw_idx=find([raw_fstr.waveform{:,2}]==wfstats(wfidx,1));
                xc.label{lblidx,3}=raw_fstr.waveform{wfidx,4};
            end
        end
    else
        continue
    end
    sums(end+1,:)={i,folder,sustIds,transIds,xc};
    save('XCORR_delay.mat','sums')
    
end
% save('x_corr_delay.mat','sufs','tfs','counter')


function [avail,out]=pre_process(folder,sustIds,transIds)
sps=30000;
trials=clearBadPerf(h5read(fullfile(folder,'events.hdf5'),'/trials')');
if isempty(trials)
    avail=false;
    out=[];
    return
end

%     trials=double(trials);
%     info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7),trials(:,8)];

% s1s=30000;
% FR_Th=1.0;

% metaf=ls(fullfile(rootpath,'*.meta'));
% fh=fopen(fullfile(rootpath,metaf));
% ts=textscan(fh,'%s','Delimiter',{'\n'});
% nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
% spkNThresh=nSample/385/s1s/2*FR_Th;
% clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
% waveformGood=strcmp(clusterInfo{:,4},'good');
% freqGood=clusterInfo{:,10}>spkNThresh;
% cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1))';
cluster_ids=[sustIds;transIds];

%  single-unit candidate


spkTS=readNPY(fullfile(folder,'spike_times.npy'));
spkId=readNPY(fullfile(folder,'spike_clusters.npy'));

FT_SPIKE=struct();

FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
for i=1:numel(cluster_ids)
    FT_SPIKE.timestamp{i}=spkTS(spkId==i)';
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



function out=clearBadPerf(facSeq)

if length(facSeq)>=40
    facSeq(:,9)=0;
    i=40;
    while i<=length(facSeq)
        goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
        if goodOff>=30 %.75 correct rate
            facSeq(i-39:i,9)=1;
        end
        i=i+1;
    end
    out=facSeq(facSeq(:,9)==1,:);
else
    out=[];
end
end


function Xc=plotxcorr(spikeTrials)
cfg             = [];
cfg.maxlag      = 0.02; % maximum 50 ms
cfg.binsize     = 0.0005; % bins of 0.5 ms
cfg.outputunit  = 'raw'; % make unit area
cfg.latency     = [1 5];
cfg.vartriallen = 'yes'; % do not allow variable trial lengths
cfg.method      = 'xcorr'; % compute the normal cross-correlogram
cfg.debias      = 'no';
Xc = ft_spike_xcorr(cfg,spikeTrials);

% compute the shuffled correlogram
% cfg.method      = 'shiftpredictor'; % compute the shift predictor
% Xshuff = ft_spike_xcorr(cfg,spikeTrials);

end

