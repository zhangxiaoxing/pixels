addpath('npy-matlab-master\npy-matlab')
addpath('fieldtrip-20200320')
ft_defaults
rootlist=dir('K:\neupix\DataSum\**\spike_times.npy');
fileCounter=1;
for f=rootlist'
    %     disp(f.folder)
    %TODO: check existence
    [avail,FT_SPIKE]=pre_process(f);
    if ~avail
        continue
    end
    %     plotPsth(FT_SPIKE);
    for i=1:20
        for j=(i+1):20
            tic
            fh=figure('color','w','Position',[100,100,2000,1200]);
            plotJPSTH(FT_SPIKE,[i,j],4,1);
            plotJPSTH(FT_SPIKE,[i,j],8,2);
            print(fh,'-dpng',sprintf('jpsth_session_%d_su_%d_%d.png',fileCounter,i,j));
            close(fh);  
            toc
            error('test unit')
        end
    end
    
end

function [avail,out]=pre_process(f)
sps=30000;
disp(f.folder)
rootpath=f.folder;
trials=clearBadPerf(h5read(fullfile(rootpath,'events.hdf5'),'/trials')');
if isempty(trials)
    avail=false;
    out=[];
    return
end

%     trials=double(trials);
%     info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7),trials(:,8)];

s1s=30000;
FR_Th=1.0;

metaf=ls(fullfile(rootpath,'*.meta'));
fh=fopen(fullfile(rootpath,metaf));
ts=textscan(fh,'%s','Delimiter',{'\n'});
nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
spkNThresh=nSample/385/s1s/2*FR_Th;
clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
waveformGood=strcmp(clusterInfo{:,4},'good');
freqGood=clusterInfo{:,10}>spkNThresh;
cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1))';
if numel(cluster_ids)<2
    avail=false;
    out=[];
    return
end
%  single-unit candidate


spkTS=readNPY(fullfile(rootpath,'spike_times.npy'));
spkId=readNPY(fullfile(rootpath,'spike_clusters.npy'));

FT_SPIKE=struct();

FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids')));
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

avail=true;
out=FT_SPIKE;

end


function plotPsth(spikeTrials)
cfg         = [];
cfg.binsize =  [0.25];
cfg.latency = [-1 4];
psth = ft_spike_psth(cfg,spikeTrials);

cfg              = [];
cfg.topplotfunc  = 'line'; % plot as a line
cfg.spikechannel = spikeTrials.label([1]);
cfg.latency      = [-1 4];
cfg.errorbars    = 'std'; % plot with the standard deviation
cfg.interactive  = 'no'; % toggle off interactive mode
figure()
ft_spike_plot_raster(cfg,spikeTrials, psth)

end


function plotJPSTH(spikeTrials,suids,samp,row)
cfg         = [];
cfg.latency = [-1 3];
cfg.timwin  = [-0.025 0.025];
cfg.fsample = 200;
[sdf] = ft_spikedensity(cfg,spikeTrials);

% compute the joint psth
cfg               = [];
cfg.trials        = spikeTrials.trialinfo(:,5)==samp;
cfg.normalization = 'no';
cfg.channelcmb    = spikeTrials.label(suids);
cfg.method        = 'jpsth';
jpsth = ft_spike_jpsth(cfg,sdf);

cfg.method = 'shiftpredictor';
jpsthShuff = ft_spike_jpsth(cfg,sdf);

% subtract the predictor
jpsthSubtr = jpsth;
jpsthSubtr.jpsth = jpsth.jpsth-jpsthShuff.shiftpredictor;

subplot(2,3,row*3-2)
ft_spike_plot_jpsth(cfg,jpsth)
ax=findall(gcf,'type','axes');
ax(3).Title.String=sprintf('S%d normalized JPSTH',row);
arrayfun(@(x) xline(ax(3),x,'--w'),[0,1]);
arrayfun(@(x) yline(ax(3),x,'--w'),[0,1]);
subplot(2,3,row*3-1)
ft_spike_plot_jpsth(cfg,jpsthShuff)
ax=findall(gcf,'type','axes');
ax(3).Title.String=sprintf('S%d shuffled JPSTH',row);
arrayfun(@(x) xline(ax(3),x,'--w'),[0,1]);
arrayfun(@(x) yline(ax(3),x,'--w'),[0,1]);
subplot(2,3,row*3)
ft_spike_plot_jpsth(cfg,jpsthSubtr)
ax=findall(gcf,'type','axes');
ax(3).Title.String=sprintf('S%d non-evoked JPSTH',row);
arrayfun(@(x) xline(ax(3),x,'--w'),[0,1]);
arrayfun(@(x) yline(ax(3),x,'--w'),[0,1]);
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