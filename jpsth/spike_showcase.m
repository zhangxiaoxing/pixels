%idces from hebb_pattern_showcase.m
% keyboard
% hidx=1
addpath('fieldtrip-20200320')
ft_defaults
load hebb_pattern.mat
load 114_sorted_file_path.mat
hebbPattern=hebbPatternS2;
load hebb_showcase_stats.mat
cmp_stats=cmp_stats_s2;
%%%%%%%%%%%%%%%%%
[~,s]=sort(cmp_stats(:,2),'descend');
cmp_stats=cmp_stats(s,:);
% cmp_stats=cmp_stats(cmp_stats(:,3)>4,:);
%%%%%%%%%%%

for hidx=hrange
    ssidx=cmp_stats(hidx,1);
    sessIdx=idivide(int32(hebbPattern{ssidx,1}(1)),int32(100000));
    delay=6;
    sample=1;
    
    if isunix
        cd('~/pixels/jpsth')
        homedir='/home/zx/neupix/wyt';
    elseif ispc
        homedir='k:\neupix\wyt';
    end
    
    folder=regexp(sorted_fpath{sessIdx},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
    currmodel='selec';
    sustIds=[];nonselIds=[];
    cocountSum=[];
    fh=figure('Color','w','Position',[100,100,1600,400]);
    for ringidx=1:3
        transIds=int32(rem(unique([hebbPattern{ssidx,ringidx}]),100000))';
        msize=numel(transIds);
        [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
        if sample==1
            cfg.trials = find(spktrial.trialinfo(:,5)==4 & spktrial.trialinfo(:,8)==delay);
        else
            cfg.trials = find(spktrial.trialinfo(:,5)==8 & spktrial.trialinfo(:,8)==delay);
        end
        samples=cell(0);
        cocount=[];
        for tidx=cfg.trials(:)'
            disp(tidx);
            ts_id=[];
            for seqid=1:msize
                tsel=spktrial.trial{seqid}==tidx;
                ts=spktrial.time{seqid}(tsel);
                ts=ts(ts>1 & ts<7);
                ts(2,:)=seqid;
                ts_id=[ts_id;ts'];
            end
            [~,s]=sort(ts_id(:,1));
            ts_id=ts_id(s,:);
            ts_id_tagged=relax_tag(ts_id,msize);
            coact_count=sum(ts_id_tagged(:,3));
            if coact_count>19
                samples{end+1}=ts_id_tagged;
                cocount(end+1,:)=[coact_count,tidx];
            end
        end
        if ~isempty(cocount)
            [cctt,sidx]=max(cocount(:,1));
            plotseq(samples{sidx},msize,cocount(sidx,2),transIds,sessIdx,ringidx)
            cocountSum=[cocountSum,cctt];
            if ringidx==3
                print(fh,sprintf('hebb_spike_seq_%d_%d_%d_%d.png',...
                    cocountSum(2),cocountSum(3),cocountSum(1),ssidx),'-dpng','-r300');
                close(fh);
            end
        else
            close(fh)
            break
        end
    end
end
%
% load(fullfile(spkFolder,'spike_info.mat'));
% suid1=spike_info{1}{1};
% suid2=spike_info{1}{2};
%
% TS1=spike_info{2}{1};
% TS2=spike_info{2}{2};
% sts=[];
% for suid=1:numel(suidr1)
%     if suidr1(suid)<10000
%         t=TS1(suid1==suidr1(suid));
%         t(:,2)=suid;
%         sts=[sts;t];
%     else
%         t=TS2(suid2==suidr1(suid));
%         t(:,2)=suid;
%         sts=[sts;t];
%     end
% end
% [ts,tsidx]=sort(sts(:,1));
% sts=sts(tsidx,:);




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

function out=tag_ring(in,msize)
out=in;
out(:,3)=0;
if msize==3
    for i=1:length(in)-3
        if (isequal(in(i:i+3,2),[1 2 3 1]') || isequal(in(i:i+3,2),[2 3 1 2]') || isequal(in(i:i+3,2),[3 1 2 3]')) ...
                && max(diff(in(i:i+3,1)))<0.01 && min(diff(in(i:i+3,1)))>0.001
            out(i:i+3,3)=1;
        end
    end
elseif msize==4
    for i=1:length(in)-4
        if (isequal(in(i:i+4,2),[1 2 3 4 1]')...
                || isequal(in(i:i+4,2),[2 3 4 1 2]')...
                || isequal(in(i:i+4,2),[3 4 1 2 3]')...
                || isequal(in(i:i+4,2),[4 1 2 3 4]'))...
                && max(diff(in(i:i+4,1)))<0.01 ...
                && min(diff(in(i:i+4,1)))>0.001
            out(i:i+3,3)=1;
        end
    end
end
end

function out=relax_tag(in,msize)
out=in;
out(:,3)=0;
skiptag=0;
for i=1:length(in)
    if i<skiptag
        continue
    end
    curr_su=in(i,2);
    targets=[(curr_su+1):(curr_su+msize-1),curr_su];
    targets(targets>msize)=targets(targets>msize)-msize;
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

function plotseq(in,msize,tidx,ids,sessIdx,ringIdx)
subplot(1,3,ringIdx)
hold on
color={'r','b','c','m'};
for bin=1:6
    for i=1:msize
        j=i;
        selp=in(:,2)==i & in(:,3)==1 & in(:,1)>=bin & in(:,1)<bin+1;
        seln=in(:,2)==i & in(:,3)==0 & in(:,1)>=bin & in(:,1)<bin+1;
        plot(repmat(in(selp,1)',2,1)-bin,(bin-1)*(msize+1)+repmat([j-0.6;j+0.6],1,nnz(selp)),'-','Color',color{i});
        plot(repmat(in(seln,1)',2,1)-bin,(bin-1)*(msize+1)+repmat([j-0.6;j+0.6],1,nnz(seln)),'-','Color',[0.8,0.8,0.8]);
    end
end
xlabel('Time within time-bin(s)')
ylabel('Time-bin')
set(gca,'XTick',0:0.5:1,'YTick',2:(msize+1):((msize+1)*6),'YTickLabel',1:6)
if msize==3
    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d',sessIdx,tidx,ids(1),ids(2),ids(3)));
    %     print(fh,sprintf('spike_seq_%d_%d_%d.png',ids(1),ids(2),ids(3)),'-dpng');
elseif msize==4
    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d, %d',sessIdx,tidx,ids(1),ids(2),ids(3),ids(4)));
    %     print(fh,sprintf('spike_seq_%d_%d_%d_%d.png',ids(1),ids(2),ids(3),ids(4)),'-dpng');
    
end
end

