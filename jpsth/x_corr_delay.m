%  A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
%  chan2. Thus, a negative lag represents a spike in the second dimension of
%  stat.xcorr before the channel in the third dimension of stat.stat.
cd('~/pixels/jpsth')
homedir='/home/zx/neupix/wyt';
currmodel='nonsel';
prefix='0618';
delay=6;
bin_range=[6 7];
addpath(fullfile('npy-matlab-master','npy-matlab'))
addpath('fieldtrip-20200320')
ft_defaults
sus_trans=h5read('../transient_6.hdf5','/sus_trans'); %export_arr = np.vstack((sust, transient, switched, unclassified, early_in_6s, late_in_6s, prefer_s))
reg_list=h5read('../transient_6.hdf5','/reg');
cid_list=h5read('../transient_6.hdf5','/cluster_id');
path_list=h5read('../transient_6.hdf5','/path');
if startsWith(currmodel,'selec')
    sust=find(sus_trans(:,1));
    trans=find(sus_trans(:,2));
    supool=[sust;trans]';
elseif startsWith(currmodel,'nonsel')
    sust=[];
    trans=[];
    nonsel_logic=~(sus_trans(:,1) | sus_trans(:,2) | sus_trans(:,3)|sus_trans(:,4));
    supool=find(nonsel_logic);
    keyboard
end
counter=[];
done=[];

error_list=cell(0);

for i=1:length(supool)
    if ismember(supool(i),done)
        continue
    end

    folder=regexp(path_list{supool(i)},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list);
    if folderType<0
        continue
    end
    wffile=fullfile(metaFolder,'wf_stats.hdf5');% posix
    if isfile(wffile)
        if startsWith(currmodel,'selec')
            if folderType==1
                sustIds=cid_list(startsWith(path_list,folder) & sus_trans(:,1));
                transIds=cid_list(startsWith(path_list,folder) & sus_trans(:,2));
                sameFolder=find(startsWith(path_list,folder) & (sus_trans(:,1)| sus_trans(:,2)));
                done=[done;sameFolder];
                sustCount=numel(sustIds);
                transCount=numel(transIds);
            elseif folderType==2
                fimec0=replace(folder,'imec1','imec0');
                fimec1=replace(folder,'imec0','imec1');

                sustIds0=cid_list(startsWith(path_list,fimec0) & sus_trans(:,1));
                transIds0=cid_list(startsWith(path_list,fimec0) & sus_trans(:,2));

                sustIds1=cid_list(startsWith(path_list,fimec1) & sus_trans(:,1))+10000;
                transIds1=cid_list(startsWith(path_list,fimec1) & sus_trans(:,2))+10000;

                sameFolder=find((startsWith(path_list,fimec0)|startsWith(path_list,fimec1)) & (sus_trans(:,1)| sus_trans(:,2)));
                done=[done;sameFolder];

                sustIds=[sustIds0(:);sustIds1(:)];
                transIds=[transIds0(:);transIds1(:)];

                sustCount=numel(sustIds);
                transCount=numel(transIds);            

            end

            if transCount<1
                continue
            end
            nonselIds=[];
            [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
        elseif startsWith(currmodel,'nonsel')
            if folderType==1
                sameFolder=find(startsWith(path_list,folder) & nonsel_logic);
                nonselIds=cid_list(startsWith(path_list,folder) & nonsel_logic);
                done=[done;sameFolder];
            elseif folderType==2
                fimec0=replace(folder,'imec1','imec0');
                fimec1=replace(folder,'imec0','imec1');

                nonselIds0=cid_list(startsWith(path_list,fimec0) & nonsel_logic);

                nonselIds1=cid_list(startsWith(path_list,fimec1) & nonsel_logic)+10000;


                sameFolder=find((startsWith(path_list,fimec0)|startsWith(path_list,fimec1)) & nonsel_logic);
                done=[done;sameFolder];

                nonselIds=[nonselIds0(:);nonselIds1(:)];
                sustIds=[];
                transIds=[];

            end

            if numel(nonselIds)<1
                continue
            end
            [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
        
        end
        if avail
            [xc_s1,xcshuf_s1,xc_s2,xcshuf_x2]=plotxcorr(spktrial,delay,bin_range);
        end
    else
        continue
    end
    sums={i,folder,sustIds,transIds,xc_s1,xcshuf_s1,xc_s2,xcshuf_x2}; %per folder save
    save(sprintf('%s_%s_XCORR_duo_f%d_delay_%d_%d_%d_2msbin.mat',prefix,currmodel,i,delay,bin_range(1),bin_range(2)),'sums','-v7.3','sust','trans','supool','counter','done') %prefix
 	fprintf('%d of %d\n',i,length(supool))
end

return 
function [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list)
    metaFolder=replace(folder,'\','/');
    metaFolder=fullfile('/home/zx/neupix/wyt/DataSum',metaFolder);
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
    s1s=30000;
    FR_Th=1.0;

    metaf=strtrim(ls(fullfile(metaFolder,'*.meta')));
    fh=fopen(metaf);
    ts=textscan(fh,'%s','Delimiter',{'\n'});
    nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
    spkNThresh=nSample/385/s1s/2*FR_Th;
    clusterInfo = readtable(fullfile(metaFolder,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
    waveformGood=strcmp(clusterInfo{:,4},'good');
    freqGood=clusterInfo{:,10}>spkNThresh;
    cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
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



function out=clearBadPerf(facSeq, model)
if strcmp(model, 'error')
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
end


function [Xc_S1,Xshuff_S1,Xc_S2,Xshuff_S2]=plotxcorr(spikeTrials,delay,bin_range)
% https://www.nature.com/articles/nn799
% A role for inhibition in shaping the temporal flow of information in prefrontal cortex
% Christos Constantinidis, Graham V. Williams & Patricia S. Goldman-Rakic 
% Nature Neuroscience volume 5, pages175-180(2002)
% 
% Neuron, Volume 76
% Functional Microcircuit Recruited during Retrieval of Object Association Memory in Monkey Perirhinal Cortex
% Toshiyuki Hirabayashi, Daigo Takeuchi, Keita Tamura, and Yasushi Miyashita


cfg             = [];
cfg.maxlag      = 0.1; % maximum 100 ms
cfg.binsize     = 0.002; % bins of 2 ms
cfg.outputunit  = 'raw'; % make unit area
cfg.latency     = bin_range; % time bin based on sample onset
cfg.vartriallen = 'no'; % allow variable trial lengths
cfg.debias      = 'no';

cfg.trials      = find(spikeTrials.trialinfo(:,5)==4 & spikeTrials.trialinfo(:,8)==delay);
if numel(cfg.trials)<2
    Xc_S1=[];
    Xshuff_S1=[];
else
    cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc_S1 = ft_spike_xcorr(cfg,spikeTrials);
    cfg.method      = 'shiftpredictor'; % compute the shift predictor
    Xshuff_S1 = ft_spike_xcorr(cfg,spikeTrials);
end

cfg.trials      = find(spikeTrials.trialinfo(:,5)==8 & spikeTrials.trialinfo(:,8)==delay);
if numel(cfg.trials)<2
    Xc_S2=[];
    Xshuff_S2=[];
else
    cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc_S2 = ft_spike_xcorr(cfg,spikeTrials);
    cfg.method      = 'shiftpredictor'; % compute the shift predictor
    Xshuff_S2 = ft_spike_xcorr(cfg,spikeTrials);
end
% compute the shuffled correlogram
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
