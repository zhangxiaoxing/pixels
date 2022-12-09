%% load
fc.util.loaded

%% dependency
% fc.util.dependency
ephys.util.dependency() %TODO determine which lib is necessary

%%
trials=behav.procPerf(h5read(fullfile(folder,'events.hdf5'),'/trials')');
fstr=load(fullfile(folder,'spike_info.mat'));
spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);

FT_SPIKE=struct();
FT_SPIKE.label=strtrim(cellstr(num2str(cids')));
FT_SPIKE.timestamp=cell(1,2);
for i=1:numel(cids)
    FT_SPIKE.timestamp{i}=spkTS(spkId==cids(i))';
end
%  continuous format F T struct file
sps=30000;
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;
FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

fh=figure('Position',[32,32,800,600],'Color','w');

subIdx=1;
for trlIdx=ptrials

    ts1=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==trlIdx)';
    ts2=FT_SPIKE.time{2}(FT_SPIKE.trial{2}==trlIdx)';
    
    ts1(:,2)=1;
    ts2(:,2)=2;
    ts_id=[ts1;ts2];
    
    [~,s]=sort(ts_id(:,1));
    ts_id=ts_id(s,:);
    ts_id_tagged=fc.fc_tag(ts_id);

    subplot(2,1,subIdx)
    hold on
    sel1bk=ts_id_tagged(:,2)==1 & ts_id_tagged(:,3)==0;
    sel2bk=ts_id_tagged(:,2)==2 & ts_id_tagged(:,3)==0;
    sel1fc=ts_id_tagged(:,2)==1 & ts_id_tagged(:,3)==1;
    sel2fc=ts_id_tagged(:,2)==2 & ts_id_tagged(:,3)==1;
    
    plot([1;1]*ts_id_tagged(sel1bk,1)',[0;1]*ones(1,nnz(sel1bk)),'k-','Color',[0.8,0.8,0.8])
    plot([1;1]*ts_id_tagged(sel2bk,1)',[1;2]*ones(1,nnz(sel2bk)),'k-','Color',[0.8,0.8,0.8])
    plot([1;1]*ts_id_tagged(sel1fc,1)',[0;1]*ones(1,nnz(sel1fc)),'r-')
    plot([1;1]*ts_id_tagged(sel2fc,1)',[1;2]*ones(1,nnz(sel2fc)),'r-')
    xlim([3,4]);
    subIdx=subIdx+1;
end

