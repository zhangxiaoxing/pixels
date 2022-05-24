if ~exist('sums_all','var')
%     load(fullfile('bzdata','sums_ring_stats_4.mat'),'sums4');
%     load(fullfile('bzdata','sums_ring_stats_3.mat'),'sums3');
    fstr=load(fullfile('bzdata','sums_ring_stats_all.mat'));    
    sums_all=[fstr.sums_all{1};fstr.sums_all{2};fstr.sums_all{3}];
end


meta=ephys.util.load_meta();
% [~,~,ratiomap]=ref.get_pv_sst();
load('OBM1Map.mat','OBM1map');
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
span=[];
starts_ends=[];
usess=unique(cell2mat(sums_all(:,1)));


for si=reshape(usess,1,[])
%     load(fullfile(fl(fi).folder,fl(fi).name),'sums');
    sums=sums_all(cell2mat(sums_all(:,1))==si,:);
    sess=sums{1,1};
    sesssel=meta.sess==sess;
    cids=meta.allcid(sesssel);
    regs=meta.reg_tree(5,sesssel);
    for ri=1:size(sums)
        ring_reg=arrayfun(@(x) regs(cids==x),sums{ri,3},'UniformOutput',false);
%         if ~all(cellfun(@(x) ratiomap.isKey(x),ring_reg),'all')
        if ~all(cellfun(@(x) OBM1map.isKey(x),ring_reg),'all')
            continue
        end
%         pv_ratio=cellfun(@(x) ratiomap(char(x)),ring_reg);
        [su_fill,ccfid_fill,hier_idx_fill]=deal(nan(1,5));
        hier_idx=cellfun(@(x) OBM1map(char(x)),ring_reg);
        ccfid=cellfun(@(x) idmap.reg2ccfid(char(x)),ring_reg);
        ccfid_fill(1:numel(ccfid))=ccfid;
        su_fill(1:numel(sums{ri,3}))=sums{ri,3};
        hier_idx_fill(1:numel(hier_idx))=hier_idx;
        span=[span;sess,su_fill,max(hier_idx)-min(hier_idx)];
        if max(hier_idx)-min(hier_idx) >0
            [SB,SBG]=groupcounts([sums{ri,5}.starts.';(1:5).']);
            [EB,EBG]=groupcounts([sums{ri,5}.ends.';(1:5).']);
            starts_ends=[starts_ends;SB.'-1,EB.'-1,ccfid_fill,hier_idx_fill];
        end
    end
end
fh=figure('Color','w');
histogram(span(:,7),0:0.5:12,'FaceColor','w','EdgeColor','k','Normalization','probability')
xlabel('OB-M1 hierarchy index span');
ylabel('Loop activity probility');

starts=struct();
[starts.HighRatio,starts.LowRatio]=deal([]);
ends=starts;
sucount=starts;
for ri=1:size(starts_ends,1)
    %% high
    maxv=nanmax(starts_ends(ri,16:20));
    hi_all=find(starts_ends(ri,16:20)==maxv);
    maxc=max(starts_ends(ri,hi_all));
    starts.HighRatio=[starts.HighRatio;maxc./sum(starts_ends(ri,1:5),'all')];
    sucount.HighRatio=[sucount.HighRatio;numel(hi_all)./nnz(isfinite(starts_ends(ri,16:20)))];
    maxc=max(starts_ends(ri,hi_all+5));
    ends.HighRatio=[ends.HighRatio;maxc./sum(starts_ends(ri,6:10),'all')];
    %% low
    minv=nanmin(starts_ends(ri,16:20));
    low_all=find(starts_ends(ri,16:20)==minv);
    maxc=max(starts_ends(ri,low_all));
    starts.LowRatio=[starts.LowRatio;maxc./sum(starts_ends(ri,1:5),'all')];
    sucount.LowRatio=[sucount.LowRatio;numel(low_all)./nnz(isfinite(starts_ends(ri,16:20)))];
    maxc=max(starts_ends(ri,low_all+5));
    ends.LowRatio=[ends.LowRatio;maxc./sum(starts_ends(ri,6:10),'all')];
end

cp=permutation_test_1d(sucount.HighRatio,sucount.LowRatio,1000);
sp=permutation_test_1d(starts.HighRatio,starts.LowRatio,10000);
ep=permutation_test_1d(ends.HighRatio,ends.LowRatio,10000);


fh=figure('Color','w','Position',[32,32,350,175]);
subplot(1,3,1)
hold on
swxx=[zeros(numel(sucount.LowRatio),1);ones(numel(sucount.HighRatio),1)]+1;
swyy=[sucount.LowRatio;sucount.HighRatio];
swarmchart(swxx(1:10:numel(swxx)),swyy(1:10:numel(swyy)),1,'o','MarkerFaceColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none','XJitterWidth',0.5)
mm=[mean(sucount.LowRatio),mean(sucount.HighRatio)]
ci=[bootci(500,@(x) mean(x),sucount.LowRatio),bootci(500,@(x) mean(x),sucount.HighRatio)]
errorbar(1:2,mm,ci(1,:)-mm,ci(2,:)-mm,'k.','CapSize',15);
set(gca(),'XTick',1:2,'XTickLabel',{'Sens.','Motor'},'XTickLabelRotation',90,'YTick',0:0.5:1,'YTickLabel',0:50:100)
ylim([0,1])
ylabel('Fraction of neurons in loops (%)');
text(1.5,1,sprintf('p=%.3f',cp),'HorizontalAlignment','center','VerticalAlignment','top');
xlim([0.5,2.5])

subplot(1,3,2)
hold on
swxx=[zeros(numel(starts.LowRatio),1);ones(numel(starts.HighRatio),1)]+1;
swyy=[starts.LowRatio;starts.HighRatio];
swarmchart(swxx(1:10:numel(swxx)),swyy(1:10:numel(swyy)),1,'o','MarkerFaceColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none','XJitterWidth',0.5)
% boxplot([starts.HighRatio;starts.LowRatio],[zeros(numel(starts.HighRatio),1);ones(numel(starts.LowRatio),1)]+1,'Colors','k','Symbol','','Widths',0.8)
mm=[mean(starts.LowRatio),mean(starts.HighRatio)]
ci=[bootci(500,@(x) mean(x),starts.LowRatio),bootci(500,@(x) mean(x),starts.HighRatio)]
errorbar(1:2,mm,ci(1,:)-mm,ci(2,:)-mm,'k.','CapSize',15);
set(gca(),'XTick',1:2,'XTickLabel',{'Sens.','Motor'},'XTickLabelRotation',90,'YTick',0:0.5:1,'YTickLabel',0:50:100)
ylim([0,1])
ylabel('Initiate a new loop activity (%)');
text(1.5,1,sprintf('p=%.3f',sp),'HorizontalAlignment','center','VerticalAlignment','top');
xlim([0.5,2.5])

subplot(1,3,3)
hold on
swxx=[zeros(numel(ends.LowRatio),1);ones(numel(ends.HighRatio),1)]+1;
swyy=[ends.LowRatio;ends.HighRatio];
swarmchart(swxx(1:10:numel(swxx)),swyy(1:10:numel(swyy)),1,'o','MarkerFaceColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none','XJitterWidth',0.5)
% boxplot([ends.HighRatio;ends.LowRatio],[zeros(numel(ends.HighRatio),1);ones(numel(ends.LowRatio),1)]+1,'Colors','k','Symbol','','Widths',0.8)
mm=[mean(ends.LowRatio),mean(ends.HighRatio)]
ci=[bootci(500,@(x) mean(x),ends.LowRatio),bootci(500,@(x) mean(x),ends.HighRatio)]
errorbar(1:2,mm,ci(1,:)-mm,ci(2,:)-mm,'k.','CapSize',15);
set(gca(),'XTick',1:2,'XTickLabel',{'Sens.','Motor'},'XTickLabelRotation',90,'YTick',0:0.5:1,'YTickLabel',0:50:100)
ylim([0,1])
ylabel('Break current loop activity (%)');
text(1.5,1,sprintf('p=%.3f',ep),'HorizontalAlignment','center','VerticalAlignment','top');
xlim([0.5,2.5])
exportgraphics(fh,'loops_starts_ends_hier.pdf','ContentType','vector')


mean([starts.HighRatio./sucount.HighRatio,starts.LowRatio./sucount.LowRatio])
mean([ends.HighRatio./sucount.HighRatio,ends.LowRatio./sucount.LowRatio])
ranksum(starts.HighRatio./sucount.HighRatio,starts.LowRatio./sucount.LowRatio)
ranksum(ends.HighRatio./sucount.HighRatio,ends.LowRatio./sucount.LowRatio)