
if ~exist('fidx','var') || ~isfile(fullfile('bzdata',sprintf('0203_BZ_XCORR_duo_f%d.mat',fidx)))
    disp('Error loading file')
    if isunix
        quit(0);
    else
        return
    end
end

trialthresh=20;

disp(fidx)
load(fullfile('bzdata',sprintf('0203_BZ_XCORR_duo_f%d.mat',fidx)),'mono','folder')

ts_sep=0;
pre_thresh=0;
trl_thresh=5;
%%
[avail,folder]=fc.util.loaded(mono,folder);
fc.util.dependency
%%

trials=h5read(fullfile(folder,'events.hdf5'),'/trials')';
% trials=behav.procPerf(trials);

sel_6s_S1=trials(:,5)==4 & trials(:,8)== 6;
sel_6s_S2=trials(:,5)==8 & trials(:,8)== 6;

sel_3s_S1=trials(:,5)==4 & trials(:,8)== 3;
sel_3s_S2=trials(:,5)==8 & trials(:,8)== 3;

if min([nnz(sel_6s_S1),nnz(sel_6s_S2),nnz(sel_3s_S1),nnz(sel_3s_S2)])<trialthresh
    disp('Performance criteria')
    quit(0);
end

sums=cell(0);
fstr=load(fullfile(folder,'spike_info.mat'));
spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);

for idx=1:size(mono.sig_con,1)
    FT_SPIKE=struct();
    cluster_ids=mono.completeIndex(mono.sig_con(idx,:),2);
    FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
    FT_SPIKE.timestamp=cell(1,2);
    for i=1:numel(cluster_ids)
        FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
    end
    %  continuous format F T struct file
    sps=30000;
    cfg=struct();
    cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
    cfg.trlunit='timestamps';
    cfg.timestampspersecond=sps;
    FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
    
    %dimord={ts1,ts2, fcevents, anti-causal fc, 50ms-lag fc} x trl(~240) x bin(-3:11)
    
    stats=zeros(5,size(trials,1),14);
    for trlIdx=1:size(trials,1)
        ts1=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==trlIdx)';
        ts2=FT_SPIKE.time{2}(FT_SPIKE.trial{2}==trlIdx)';
        ts1=ts1([false;diff(ts1)>ts_sep]);
        if isempty(ts1)
            continue
        end
        ts1(:,2)=1;
        ts2(:,2)=2;
        ts_id=[ts1;ts2];
        [~,s]=sort(ts_id(:,1));
        ts_id=ts_id(s,:);
        ts_id_tagged=fc.fc_tag(ts_id,false);
        stats(1,trlIdx,:)=histcounts(ts1(:,1),-3:11);
        stats(2,trlIdx,:)=histcounts(ts2(:,1),-3:11);
        stats(3,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,3)==1 & ts_id_tagged(:,2)==1,1),-3:11);
        stats(4,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,4)==1 & ts_id_tagged(:,2)==1,1),-3:11);
        stats(5,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,5)==1 & ts_id_tagged(:,2)==1,1),-3:11);
    end
    sums(end+1,:)={idx,cluster_ids,stats};
    if exist('debug','var') && debug && size(sums,1)>=50
        break
    end
end
save(sprintf('fc_decoding_f%d.mat',fidx),'sums','trials','folder')
if isunix
    quit(0)
else
    return
end



