fl=dir('bzdata\ring_stats_3_*.mat');
meta=ephys.util.load_meta();
[~,~,ratiomap]=ref.get_pv_sst();
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
span=[];
starts_ends=[];
for fi=1:size(fl,1)
    load(fullfile(fl(fi).folder,fl(fi).name),'sums');
    sess=sums{1,1};
    sesssel=meta.sess==sess;
    cids=meta.allcid(sesssel);
    regs=meta.reg_tree(5,sesssel);
    for ri=1:size(sums)
        ring_reg=arrayfun(@(x) regs(cids==x),sums{ri,3},'UniformOutput',false);
        if ~all(cellfun(@(x) ratiomap.isKey(x),ring_reg),'all')
            continue
        end
        pv_ratio=cellfun(@(x) ratiomap(char(x)),ring_reg);
        ccfid=cellfun(@(x) idmap.reg2ccfid(char(x)),ring_reg);
        span=[span;sess,sums{ri,3},max(pv_ratio)-min(pv_ratio)];
        if max(pv_ratio)-min(pv_ratio) >0
            [SB,SBG]=groupcounts(sums{ri,5}.starts.');
            [EB,EBG]=groupcounts(sums{ri,5}.ends.');
            if ~isequal(SBG.',1:3) || ~isequal(EBG.',1:3)
                keyboard()
            end
            starts_ends=[starts_ends;SB.',EB.',ccfid,pv_ratio];
        end
    end
end
fh=figure('Color','w');
histogram(span(:,5),0:0.05:0.5,'FaceColor','w','EdgeColor','k','Normalization','probability')
xlabel('PV-SST ratio span');
ylabel('Recycle activity probility');

% [starts.HighMost,starts.HighLeast,starts.LowMost,starts.LowLeast]=deal(0);
% ends=starts;
% for ri=1:size(starts_ends,1)
%     [~,low]=max(starts_ends(ri,7:9));
%     [~,high]=min(starts_ends(ri,7:9));
%     [~,sMost]=max(starts_ends(ri,1:3));
%     [~,eMost]=max(starts_ends(ri,4:6));
%     [~,sLeast]=min(starts_ends(ri,1:3));
%     [~,eLeast]=min(starts_ends(ri,4:6));
%     
%     if low==sMost, starts.LowMost=starts.LowMost+1; end
%     if low==sLeast, starts.LowLeast=starts.LowMost+1; end
%     if high==sMost, starts.HighMost=starts.HighMost+1; end
%     if high==sLeast, starts.HighLeast=starts.HighLeast+1; end
% 
%     if low==eMost, ends.LowMost=ends.LowMost+1; end
%     if low==eLeast, ends.LowLeast=ends.LowMost+1; end
%     if high==eMost, ends.HighMost=ends.HighMost+1; end
%     if high==eLeast, ends.HighLeast=ends.HighLeast+1; end
% end

starts=struct();
[starts.HighRatio,starts.LowRatio]=deal([]);
ends=starts;
for ri=1:size(starts_ends,1)
    if numel(unique(starts_ends(ri,7:9)))==3 % 3 regions
        [~,hh]=max(starts_ends(ri,10:12));
        starts.HighRatio=[starts.HighRatio;starts_ends(ri,hh)./sum(starts_ends(ri,1:3),'all')];
        ends.HighRatio=[ends.HighRatio;starts_ends(ri,hh+3)./sum(starts_ends(ri,4:6),'all')];
        [~,ll]=min(starts_ends(ri,10:12));
        starts.LowRatio=[starts.LowRatio;starts_ends(ri,ll)./sum(starts_ends(ri,1:3),'all')];
        ends.LowRatio=[ends.LowRatio;starts_ends(ri,ll+3)./sum(starts_ends(ri,4:6),'all')];
    else %2 regions
        if diff(maxk(starts_ends(ri,10:12),2))==0 % 2xhigh
            [~,hh]=maxk(starts_ends(ri,10:12),2);
            if starts_ends(ri,hh(1))>starts_ends(ri,hh(2))
                starts.HighRatio=[starts.HighRatio;starts_ends(ri,hh(1))./sum(starts_ends(ri,1:3),'all')];
            else
                starts.HighRatio=[starts.HighRatio;starts_ends(ri,hh(2))./sum(starts_ends(ri,1:3),'all')];
            end
            if starts_ends(ri,hh(1)+3)>starts_ends(ri,hh(2)+3)
                ends.HighRatio=[ends.HighRatio;starts_ends(ri,hh(1)+3)./sum(starts_ends(ri,4:6),'all')];
            else
                ends.HighRatio=[ends.HighRatio;starts_ends(ri,hh(2)+3)./sum(starts_ends(ri,4:6),'all')];
            end
            
            [~,ll]=min(starts_ends(ri,10:12));
            starts.LowRatio=[starts.LowRatio;starts_ends(ri,ll)./sum(starts_ends(ri,1:3),'all')];
            ends.LowRatio=[ends.LowRatio;starts_ends(ri,ll+3)./sum(starts_ends(ri,4:6),'all')];
        else %2xlow
            [~,hh]=max(starts_ends(ri,10:12));
            starts.HighRatio=[starts.HighRatio;starts_ends(ri,hh)./sum(starts_ends(ri,1:3),'all')];
            ends.HighRatio=[ends.HighRatio;starts_ends(ri,hh+3)./sum(starts_ends(ri,4:6),'all')];
            [~,ll]=mink(starts_ends(ri,10:12),2);
            if starts_ends(ri,ll(1))>starts_ends(ri,ll(2))
                starts.LowRatio=[starts.LowRatio;starts_ends(ri,ll(1))./sum(starts_ends(ri,1:3),'all')];
            else
                starts.LowRatio=[starts.LowRatio;starts_ends(ri,ll(2))./sum(starts_ends(ri,1:3),'all')];
            end
            if starts_ends(ri,ll(1)+3)>starts_ends(ri,ll(2)+3)
                ends.LowRatio=[ends.LowRatio;starts_ends(ri,ll(1)+3)./sum(starts_ends(ri,4:6),'all')];
            else
                ends.LowRatio=[ends.LowRatio;starts_ends(ri,ll(2)+3)./sum(starts_ends(ri,4:6),'all')];
            end
        end
    end
end

sp=ranksum(starts.HighRatio,starts.LowRatio);
ep=ranksum(ends.HighRatio,ends.LowRatio);

fh=figure('Color','w');
subplot(1,2,1)
hold on
swarmchart([zeros(numel(starts.HighRatio),1);ones(numel(starts.LowRatio),1)]+1,[starts.HighRatio;starts.LowRatio],9,'o','MarkerFaceColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeColor','none')
boxplot([starts.HighRatio;starts.LowRatio],[zeros(numel(starts.HighRatio),1);ones(numel(starts.LowRatio),1)]+1,'Colors','k','Symbol','','Widths',0.8)
set(gca(),'XTick',1:2,'XTickLabel',{'Higher','Lower'},'XTickLabelRotation',90,'YTick',0:0.5:1,'YTickLabel',0:50:100)
ylim([0,1])
ylabel('Probablity of starting a new recurrent activity');
text(1.5,1,sprintf('p=%.3f',sp),'HorizontalAlignment','center','VerticalAlignment','top');
subplot(1,2,2)
hold on
swarmchart([zeros(numel(ends.HighRatio),1);ones(numel(ends.LowRatio),1)]+1,[ends.HighRatio;ends.LowRatio],9,'o','MarkerFaceColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeColor','none')
boxplot([ends.HighRatio;ends.LowRatio],[zeros(numel(ends.HighRatio),1);ones(numel(ends.LowRatio),1)]+1,'Colors','k','Symbol','','Widths',0.8)
set(gca(),'XTick',1:2,'XTickLabel',{'Higher','Lower'},'XTickLabelRotation',90,'YTick',0:0.5:1,'YTickLabel',0:50:100)
ylim([0,1])
ylabel('Probablity of stop current recurrent activity');
text(1.5,1,sprintf('p=%.3f',ep),'HorizontalAlignment','center','VerticalAlignment','top');