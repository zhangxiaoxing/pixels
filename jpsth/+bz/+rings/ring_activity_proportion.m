% load pstats3, pstats4 from file or function
if ~exist('pstats3','var') || ~exist('pstats4','var')
    if true
        load('loops_proportion_stats_3.mat','pstats3');
        load('loops_proportion_stats_4.mat','pstats4');
    else
        pstats3=load_data(3);
        pstats4=load_data(4);
    end
end

sess3=cellfun(@(x) str2double(replace(x,'s','')), fieldnames(pstats3));
sess4=cellfun(@(x) str2double(replace(x,'s','')), fieldnames(pstats4));
usess=unique([sess3;sess4]);

stats=[];
for sessid=usess.'
    sf=sprintf('s%d',sessid);
    [spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessid);
    if isfield(pstats3,sf)
        ucid3=cellfun(@(x) str2double(replace(x,'c','')),fieldnames(pstats3.(sf)));
    end
    if isfield(pstats4,sf)
        ucid4=cellfun(@(x) str2double(replace(x,'c','')),fieldnames(pstats4.(sf)));
    end
    ucid=unique([ucid3;ucid4]);
    
    for onecid=ucid.'
        disp([sessid,onecid]);
        cf=sprintf('c%d',onecid);
        totalspk=nnz(spkID==onecid);
        max_prop=[];
        loop_spks=[];
        if isfield(pstats3.(sf),cf)
            max_prop=[max_prop;cellfun(@(x) numel(pstats3.(sf).(cf).(x)),fieldnames(pstats3.(sf).(cf)))];
            loop_spks=[loop_spks;cell2mat(cellfun(@(x) pstats3.(sf).(cf).(x),fieldnames(pstats3.(sf).(cf)),'UniformOutput',false))];
        end
        if isfield(pstats4.(sf),cf)
            max_prop=[max_prop;cellfun(@(x) numel(pstats4.(sf).(cf).(x)),fieldnames(pstats4.(sf).(cf)))];
            loop_spks=[loop_spks;cell2mat(cellfun(@(x) pstats4.(sf).(cf).(x),fieldnames(pstats4.(sf).(cf)),'UniformOutput',false))];
        end
        stats=[stats;sessid,onecid,totalspk,numel(max_prop),max(max_prop),numel(unique(loop_spks))];
    end
end

stats(:,7)=stats(:,5)./stats(:,3).*100;
stats(:,8)=stats(:,6)./stats(:,3).*100;
save('loops_proportion_stats.mat','stats')

singleci=bootci(500,@(x) mean(x),stats(:,7));
alterci=bootci(500,@(x) mean(x),stats(:,8));
mm=mean(stats(:,7:8));

fh=figure('Color','w','Position',[32,32,195,235]);
hold on;
swarmchart([ones(size(stats,1),1);2*ones(size(stats,1),1)],[stats(:,7);stats(:,8)],1,'o','filled','MarkerFaceColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none','XJitterWidth',0.5)
% bar(1:2,mm,0.6,'FaceColor','none','EdgeColor','k','LineWidth',1);
errorbar(1:2,mm,[singleci(1),alterci(1)]-mm,[singleci(2),alterci(2)]-mm,'k.','CapSize',15);
ylabel('Loops-associated spikes (%)');
set(gca(),'XTick',1:2,'XTickLabel',{'Single loop','Alternative paths'},'XTickLabelRotation',90)
xlim([0.5,2.5]);
ylim([0,60]);
exportgraphics(fh,'loops_spk_proportion.pdf');

function pstats=load_data(rsize)
if rsize==4
    load(fullfile('bzdata','sums_ring_stats_4.mat'),'sums4');
    rstats=sums4;
else
    load(fullfile('bzdata','sums_ring_stats_3.mat'),'sums3');
    rstats=sums3;
end
usess=unique(cell2mat(rstats(:,1)));
pstats=struct();
for sessid=usess.'
    pstats.(sprintf('s%d',sessid))=struct();
    [spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessid);
    rid=find(cell2mat(rstats(:,1))==sessid);
    for ri=reshape(rid,1,[])
        disp([sessid,ri]);
        ts_id=[];
        cids=rstats{ri,3};
        per_cid_spk_cnt=cids;
        for in_ring_pos=1:numel(cids) % TODO, 1:rsize
            one_ring_sel=spkID==cids(in_ring_pos);
            per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
            ts_id=cat(1,ts_id,[spkTS(one_ring_sel),ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
        end
        ts_id=sortrows(ts_id,1);
        ts_id=[ts_id,rstats{ri,5}.tags];
        %% loops PSTH
%         ts_ring=ts_id(ts_id(:,3)~=0, 1:2);
%         SP=splitapply(@(x) {x}, ts_ring(:,1), ts_ring(:,2));
%         FT_LOOPS=struct();
%         FT_LOOPS.label=arrayfun(@(x) num2str(x),cids,'UniformOutput',false);
%         FT_LOOPS.timestamp=SP;
%         sps=30000;
%         cfg=struct();
%         cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
%         cfg.trlunit='timestamps';
%         cfg.timestampspersecond=sps;
%         ephys.util.dependency('ft',true,'buz',false); %data path and lib path dependency
%         FT_LOOPS=ft_spike_maketrials(cfg,FT_LOOPS);
%         
        %%
        for cidx=1:numel(rstats{ri,3})
            cid=rstats{ri,3}(cidx);
            if ~isfield(pstats.(sprintf('s%d',sessid)),sprintf('c%d',cid))
                pstats.(sprintf('s%d',sessid)).(sprintf('c%d',cid))=struct();
            end
            pstats.(sprintf('s%d',sessid)).(sprintf('c%d',cid)).(sprintf('r%d',ri))=ts_id(ts_id(:,2)==cidx & ts_id(:,3)~=0,1);
        end
    end
end
end