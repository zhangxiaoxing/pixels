%  A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
%  chan2. Thus, a negative lag represents a spike in the second dimension of
%  stat.xcorr before the channel in the third dimension of stat.stat.

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



function misc()
% these are helper functions that loosely related to the above functions

% for i=1:size(stats,1)
%     pci=stats{i,8};
%     if stats{i,7}<min(pci) || stats{i,7}>max(pci)
%         stats{i,11}=1;
%     else
%         stats{i,11}=0;
%     end
% end
% goodstats=stats([stats{:,11}]==1,:);
% 
% for i=1:size(goodstats,1)
%     goodstats{i,12}=sums{goodstats{i,1},2};
% end
% 
% 
% for i=1:size(goodstats,1)
%     xc=sums{goodstats{i,1},5};
%     goodstats{i,13}=xc.label{goodstats{i,2},2}(1);
%     goodstats{i,14}=xc.label{goodstats{i,3},2}(1);
% end
% 
% 
% sufs=importdata('su_list.csv',',');
% tfs=importdata('transient_6.csv');
% 
% for i=1:size(goodstats,1)
%     path=goodstats{i,12};
%     id1=goodstats{i,13};
%     id2=goodstats{i,14};
%     selid1=find(strcmp(sufs.textdata,path) & sufs.data==id1);
%     selid2=find(strcmp(sufs.textdata,path) & sufs.data==id2);
%     if ismember(1,tfs.data(8:13,selid1)) && ~ismember(2,tfs.data(8:13,selid1))
%         pref1=1;
%     elseif ismember(2,tfs.data(8:13,selid1)) && ~ismember(1,tfs.data(8:13,selid1))
%         pref1=2;
%     else
%         fprintf('%d prefered sample ambiguous')
%     end
%     
%     if ismember(1,tfs.data(8:13,selid2)) && ~ismember(2,tfs.data(8:13,selid2))
%         pref2=1;
%     elseif ismember(2,tfs.data(8:13,selid2)) && ~ismember(1,tfs.data(8:13,selid2))
%         pref2=2;
%     else
%         fprintf('%d prefered sample ambiguous')
%     end
%     goodstats{i,15}=pref1;
%     goodstats{i,16}=pref2;
%     
% end
% 
% 
% % hist, waveform
% 
% for i=1:size(goodstats,1)
%     goodstats{i,17}=squeeze(sums{goodstats{i,1},5}.xcorr(goodstats{i,2},goodstats{i,3},:));
% end
% 
% 
% for i=1:size(goodstats,1)
%     goodstats{i,18}=squeeze(sums{goodstats{i,1},5}.label{goodstats{i,2},3});
%     goodstats{i,19}=squeeze(sums{goodstats{i,1},5}.label{goodstats{i,3},3});
% end
% 
% xcorr_share=goodstats(:,[12:16,9,10,18,19,4:8,17]);
% 
% xcorr_share([xcorr_share{:,2}]<[xcorr_share{:,3}],:)=[]
% 
% save('xcorr_share_baseline.mat','xcorr_share')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for all stats



for i=1:size(y,1)
    y{i,12}=sums{y{i,1},2};
end


for i=1:size(y,1)
    xc=sums{y{i,1},5};
    y{i,13}=xc.label{y{i,2},2}(1);
    y{i,14}=xc.label{y{i,3},2}(1);
end


sufs=importdata('su_list.csv',',');
tfs=importdata('transient_6.csv');

for i=1:size(y,1)
    path=y{i,12};
    id1=y{i,13};
    id2=y{i,14};
    pref1=nan;
    pref2=nan;
    selid1=find(strcmp(sufs.textdata,path) & sufs.data==id1);
    selid2=find(strcmp(sufs.textdata,path) & sufs.data==id2);
    if ismember(1,tfs.data(8:13,selid1)) && ~ismember(2,tfs.data(8:13,selid1))
        pref1=1;
    elseif ismember(2,tfs.data(8:13,selid1)) && ~ismember(1,tfs.data(8:13,selid1))
        pref1=2;
    else
        fprintf('%d prefered sample ambiguous')
    end
    
    if ismember(1,tfs.data(8:13,selid2)) && ~ismember(2,tfs.data(8:13,selid2))
        pref2=1;
    elseif ismember(2,tfs.data(8:13,selid2)) && ~ismember(1,tfs.data(8:13,selid2))
        pref2=2;
    else
        fprintf('%d prefered sample ambiguous')
    end
    y{i,15}=pref1;
    y{i,16}=pref2;
    
end


% hist, waveform

for i=1:size(y,1)
    y{i,17}=squeeze(sums{y{i,1},5}.xcorr(y{i,2},y{i,3},:));
end


for i=1:size(y,1)
    y{i,18}=squeeze(sums{y{i,1},5}.label{y{i,2},3});
    y{i,19}=squeeze(sums{y{i,1},5}.label{y{i,3},3});
end


end





function stats()
%these are fast and dirty statistics for previous data

%xcorrstats
count_samePref_PyrPyr=sum(([x{:,4}]==[x{:,5}])' & strcmp(x(:,6),'Pyr') & strcmp(x(:,7),'Pyr'))
count_samePref_IntInt=sum(([x{:,4}]==[x{:,5}])' & strcmp(x(:,6),'Int') & strcmp(x(:,7),'Int'))
count_samePref_PyrInt=sum(([x{:,4}]==[x{:,5}])' & ~strcmp(x(:,6),x(:,7)))

count_oppoPref_PyrPyr=sum(([x{:,4}]~=[x{:,5}])' & strcmp(x(:,6),'Pyr') & strcmp(x(:,7),'Pyr'))
count_oppoPref_IntInt=sum(([x{:,4}]~=[x{:,5}])' & strcmp(x(:,6),'Int') & strcmp(x(:,7),'Int'))
count_oppoPref_PyrInt=sum(([x{:,4}]~=[x{:,5}])' & ~strcmp(x(:,6),x(:,7)))



count_samePref=sum(([x{:,4}]==[x{:,5}])')
count_oppoPref=sum(([x{:,4}]~=[x{:,5}])')

% sum(([x{:,4}]==[x{:,5}])' & strcmp(x(:,6),'Pyr') & strcmp(x(:,7),'Pyr'))
% sum(([x{:,4}]==[x{:,5}])' & strcmp(x(:,6),'Int') & strcmp(x(:,7),'Int'))
% sum(([x{:,4}]==[x{:,5}])' & ~strcmp(x(:,6),x(:,7)))
% 
% sum(([x{:,4}]~=[x{:,5}])' & strcmp(x(:,6),'Pyr') & strcmp(x(:,7),'Pyr'))
% sum(([x{:,4}]~=[x{:,5}])' & strcmp(x(:,6),'Int') & strcmp(x(:,7),'Int'))
% sum(([x{:,4}]~=[x{:,5}])' & ~strcmp(x(:,6),x(:,7)))



sum(([x{:,4}]~=[x{:,5}])' & strcmp(x(:,6),'Int') & strcmp(x(:,7),'Pyr') & [x{:,13}]'<0.5)
sum(([x{:,4}]~=[x{:,5}])' & strcmp(x(:,6),'Pyr') & strcmp(x(:,7),'Int') & [x{:,13}]'>0.5)

sum(([x{:,4}]~=[x{:,5}])' & strcmp(x(:,6),'Int') & strcmp(x(:,7),'Pyr') & [x{:,13}]'>0.5)
sum(([x{:,4}]~=[x{:,5}])' & strcmp(x(:,6),'Pyr') & strcmp(x(:,7),'Int') & [x{:,13}]'<0.5)



% all potention pairs
tot_samePref_PyrPyr=sum(([y{:,15}]==[y{:,16}])' & strcmp(y(:,9),'Pyr') & strcmp(y(:,10),'Pyr') &~strcmp(y(:,5),'auto-corr'))
tot_samePref_IntInt=sum(([y{:,15}]==[y{:,16}])' & strcmp(y(:,9),'Int') & strcmp(y(:,10),'Int') &~strcmp(y(:,5),'auto-corr'))
tot_samePref_PyrInt=sum(([y{:,15}]==[y{:,16}])' & ~strcmp(y(:,9),y(:,10)))


tot_oppoPref_PyrPyr=sum(([y{:,15}]~=[y{:,16}])' & strcmp(y(:,9),'Pyr') & strcmp(y(:,10),'Pyr') &~strcmp(y(:,5),'auto-corr'))
tot_oppoPref_IntInt=sum(([y{:,15}]~=[y{:,16}])' & strcmp(y(:,9),'Int') & strcmp(y(:,10),'Int') &~strcmp(y(:,5),'auto-corr'))
tot_oppoPref_PyrInt=sum(([y{:,15}]~=[y{:,16}])' & ~strcmp(y(:,9),y(:,10)))


tot_samePref=sum(([y{:,15}]==[y{:,16}])')
tot_oppoPref=sum(([y{:,15}]~=[y{:,16}])')


count_samePref_PyrPyr/tot_samePref_PyrPyr
count_samePref_IntInt/tot_samePref_IntInt
count_samePref_PyrInt/tot_samePref_PyrInt

count_oppoPref_PyrPyr/tot_oppoPref_PyrPyr
count_oppoPref_IntInt/tot_oppoPref_IntInt
count_oppoPref_PyrInt/tot_oppoPref_PyrInt


binofit(count_oppoPref_IntInt,tot_oppoPref_IntInt)


[tbl,chi2,pIntInt]=crosstab([zeros(tot_samePref_IntInt,1);ones(tot_oppoPref_IntInt,1)],...
    [zeros(tot_samePref_IntInt-count_samePref_IntInt,1);ones(count_samePref_IntInt,1);...
    zeros(tot_oppoPref_IntInt-count_oppoPref_IntInt,1);ones(count_oppoPref_IntInt,1)])


[tbl,chi2,pPyrPyr]=crosstab([zeros(tot_samePref_PyrPyr,1);ones(tot_oppoPref_PyrPyr,1)],...
    [zeros(tot_samePref_PyrPyr-count_samePref_PyrPyr,1);ones(count_samePref_PyrPyr,1);...
    zeros(tot_oppoPref_PyrPyr-count_oppoPref_PyrPyr,1);ones(count_oppoPref_PyrPyr,1)])

[tbl,chi2,pPyrInt]=crosstab([zeros(tot_samePref_PyrInt,1);ones(tot_oppoPref_PyrInt,1)],...
    [zeros(tot_samePref_PyrInt-count_samePref_PyrInt,1);ones(count_samePref_PyrInt,1);...
    zeros(tot_oppoPref_PyrInt-count_oppoPref_PyrInt,1);ones(count_oppoPref_PyrInt,1)])


[tbl,chi2,pSum]=crosstab([zeros(tot_samePref,1);ones(tot_oppoPref,1)],...
    [zeros(tot_samePref-count_samePref,1);ones(count_samePref,1);...
    zeros(tot_oppoPref-count_oppoPref,1);ones(count_oppoPref,1)])

[~,pciSPyrPyr]=binofit(count_samePref_PyrPyr,tot_samePref_PyrPyr)
[~,pciSPyrInt]=binofit(count_samePref_PyrInt,tot_samePref_PyrInt)
[~,pciSIntInt]=binofit(count_samePref_IntInt,tot_samePref_IntInt)
[~,pciSSum]=binofit(count_samePref,tot_samePref)

[~,pciOPyrPyr]=binofit(count_oppoPref_PyrPyr,tot_oppoPref_PyrPyr)
[~,pciOPyrInt]=binofit(count_oppoPref_PyrInt,tot_oppoPref_PyrInt)
[~,pciOIntInt]=binofit(count_oppoPref_IntInt,tot_oppoPref_IntInt)
[~,pciOSum]=binofit(count_oppoPref,tot_oppoPref)

figure('Color','w')
hold on
bar(1,count_samePref_PyrPyr/tot_samePref_PyrPyr,'w')
bar(2,count_oppoPref_PyrPyr/tot_oppoPref_PyrPyr,'k')

bar(4,count_samePref_PyrInt/tot_samePref_PyrInt,'w')
bar(5,count_oppoPref_PyrInt/tot_oppoPref_PyrInt,'k')

bar(7,count_samePref_IntInt/tot_samePref_IntInt,'w')
bar(8,count_oppoPref_IntInt/tot_oppoPref_IntInt,'k')

bhs=bar(10,count_samePref/tot_samePref,'w')
bho=bar(11,count_oppoPref/tot_oppoPref,'k')

set(gca(),'XTick',[1 4 7 10]+0.5,'XTickLabel',{'Pyr-Pyr','Pyr-Int','Int-Int','Sum'})
ylabel('functional connectivity ratio')

errorbar([1 4 7 10],[count_samePref_PyrPyr/tot_samePref_PyrPyr,...
    count_samePref_PyrInt/tot_samePref_PyrInt,...
    count_samePref_IntInt/tot_samePref_IntInt,...
    count_samePref/tot_samePref],...
    [pciSPyrPyr(2)-count_samePref_PyrPyr/tot_samePref_PyrPyr,...
    pciSPyrInt(2)-count_samePref_PyrInt/tot_samePref_PyrInt,...
    pciSIntInt(2)-count_samePref_IntInt/tot_samePref_IntInt,...
    pciSSum(2)-count_samePref/tot_samePref],'k.')

errorbar([1 4 7 10]+1,[count_oppoPref_PyrPyr/tot_oppoPref_PyrPyr,...
    count_oppoPref_PyrInt/tot_oppoPref_PyrInt,...
    count_oppoPref_IntInt/tot_oppoPref_IntInt,...
    count_oppoPref/tot_oppoPref],...
    [pciOPyrPyr(2)-count_oppoPref_PyrPyr/tot_oppoPref_PyrPyr,...
    pciOPyrInt(2)-count_oppoPref_PyrInt/tot_oppoPref_PyrInt,...
    pciOIntInt(2)-count_oppoPref_IntInt/tot_oppoPref_IntInt,...
    pciOSum(2)-count_oppoPref/tot_oppoPref],'k.')
    
legend([bhs,bho],{'same preferred','opposite preffered'})

text(1.5,0.03,sprintf('p=%0.3f',pPyrPyr),'HorizontalAlignment','center')
text(4.5,0.03,sprintf('p=%0.3f',pPyrInt),'HorizontalAlignment','center')
text(7.5,0.03,sprintf('p=%0.3f',pIntInt),'HorizontalAlignment','center')
text(10.5,0.03,sprintf('p=%0.3f',pSum),'HorizontalAlignment','center')
ylim([0,0.035])



count_S=nnz([xcorr_share{:,4}]==[xcorr_share{:,5}])
count_O=nnz([xcorr_share{:,4}]~=[xcorr_share{:,5}])

[phatS,pciS]=binofit(count_S,length(xcorr_share))
[phatO,pciO]=binofit(count_O,length(xcorr_share))

[tbl,chi2,p]=crosstab([zeros(length(xcorr_share),1);ones(length(xcorr_share),1)],...
    [zeros(length(xcorr_share)-count_S,1);ones(count_S,1);...
    zeros(length(xcorr_share)-count_O,1);ones(count_O,1)])


figure('Color','w')
hold on
bhs=bar(1,count_S/length(xcorr_share),'w');
bho=bar(2,count_O/length(xcorr_share),'k');

errorbar([1 2],[count_S/length(xcorr_share),...
    count_O/length(xcorr_share)],...
    [pciS(2)-count_S/length(xcorr_share),...
    pciO(2)-count_O/length(xcorr_share)],'k.')

legend([bhs,bho],{'same prefered sample','opposite prefered sample'})

text(1.5,0.6,'***','HorizontalAlignment','center')
ylabel('ratio of functional connectivity')
ylim([0,0.7])
set(gca(),'XTick',[])


sust_S=nnz([xcorr_share{:,4}]'==[xcorr_share{:,5}]' & (strcmp(xcorr_share(:,10),'sust') | strcmp(xcorr_share(:,11),'sust')))
sust_O=nnz([xcorr_share{:,4}]'~=[xcorr_share{:,5}]' & (strcmp(xcorr_share(:,10),'sust') | strcmp(xcorr_share(:,11),'sust')))


[phatS,pciS]=binofit(sust_S,length(xcorr_share))
[phatO,pciO]=binofit(sust_O,length(xcorr_share))

[tbl,chi2,p]=crosstab([zeros(length(xcorr_share),1);ones(length(xcorr_share),1)],...
    [zeros(length(xcorr_share)-sust_S,1);ones(sust_S,1);...
    zeros(length(xcorr_share)-sust_O,1);ones(sust_O,1)])



figure('Color','w')
hold on
bhs=bar(1,sust_S/length(xcorr_share),'k');
bho=bar(2,sust_O/length(xcorr_share),'w');

errorbar([1 2],[sust_S/length(xcorr_share),...
    sust_O/length(xcorr_share)],...
    [pciS(2)-sust_S/length(xcorr_share),...
    pciO(2)-sust_O/length(xcorr_share)],'k.')

legend([bhs,bho],{'same prefered sample','opposite prefered sample'})

text(1.5,0.025,'***','HorizontalAlignment','center')
ylabel('ratio of functional connectivity')
xlabel('sustained')
ylim([0,0.03])
set(gca(),'XTick',[])




trans_S=nnz([xcorr_share{:,4}]'==[xcorr_share{:,5}]' & (strcmp(xcorr_share(:,10),'transient') | strcmp(xcorr_share(:,11),'transient')))
trans_O=nnz([xcorr_share{:,4}]'~=[xcorr_share{:,5}]' & (strcmp(xcorr_share(:,10),'transient') | strcmp(xcorr_share(:,11),'transient')))



[phatS,pciS]=binofit(trans_S,length(xcorr_share))
[phatO,pciO]=binofit(trans_O,length(xcorr_share))

[tbl,chi2,p]=crosstab([zeros(length(xcorr_share),1);ones(length(xcorr_share),1)],...
    [zeros(length(xcorr_share)-trans_S,1);ones(trans_S,1);...
    zeros(length(xcorr_share)-trans_O,1);ones(trans_O,1)])



figure('Color','w')
hold on
bhs=bar(1,trans_S/length(xcorr_share),'k');
bho=bar(2,trans_O/length(xcorr_share),'w');

errorbar([1 2],[trans_S/length(xcorr_share),...
    trans_O/length(xcorr_share)],...
    [pciS(2)-trans_S/length(xcorr_share),...
    pciO(2)-trans_O/length(xcorr_share)],'k.')

legend([bhs,bho],{'same prefered sample','opposite prefered sample'})

text(1.5,0.6,'***','HorizontalAlignment','center')
ylabel('ratio of functional connectivity')
xlabel('transient')
ylim([0,0.7])
set(gca(),'XTick',[])






sust_S=nnz([xcorr_share{:,4}]'==[xcorr_share{:,5}]' & (strcmp(xcorr_share(:,10),'sust') & strcmp(xcorr_share(:,11),'sust')))
sust_O=nnz([xcorr_share{:,4}]'~=[xcorr_share{:,5}]' & (strcmp(xcorr_share(:,10),'sust') & strcmp(xcorr_share(:,11),'sust')))


[phatS,pciS]=binofit(sust_S,length(xcorr_share))
[phatO,pciO]=binofit(sust_O,length(xcorr_share))

[tbl,chi2,p]=crosstab([zeros(length(xcorr_share),1);ones(length(xcorr_share),1)],...
    [zeros(length(xcorr_share)-sust_S,1);ones(sust_S,1);...
    zeros(length(xcorr_share)-sust_O,1);ones(sust_O,1)])



figure('Color','w')
hold on
bhs=bar(1,sust_S/length(xcorr_share),'k');
bho=bar(2,sust_O/length(xcorr_share),'w');

errorbar([1 2],[sust_S/length(xcorr_share),...
    sust_O/length(xcorr_share)],...
    [pciS(2)-sust_S/length(xcorr_share),...
    pciO(2)-sust_O/length(xcorr_share)],'k.')

legend([bhs,bho],{'same prefered sample','opposite prefered sample'})

text(1.5,0.025,'***','HorizontalAlignment','center')
ylabel('ratio of functional connectivity')
xlabel('sustained')
ylim([0,0.03])
set(gca(),'XTick',[])


end