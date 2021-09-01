% {FWD/RWD,H2L/L2H/Local, S1/S2*correct/error*3s/6s}
% congru, incongru, nonmem
function plot_coding
% [~,com_meta]=wave.per_region_COM();
[~,~,ratiomap]=ref.get_pv_sst();
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
[metas,stats]=bz.fccoding.get_fc_coding();
% commap=containers.Map('KeyType','int32','ValueType','any');
% for ci=1:numel(com_meta,1)
%     commap(bitshift(com_meta{ci,1},16)+com_meta{ci,2})=com_meta{ci,3};
% end

% ccfids=unique(metas(:,6:7));
% for ri=1:numel(ccfids)
%     if ~ratiomap.isKey(idmap.ccfid2reg(ccfids(ri)))
%         disp(ccfids(ri))
%     end
% end

% pctmm=mean([stats(congrus1,5);stats(congrus2,6)]);
% pctci=bootci(1000,@(x) mean(x), [stats(congrus1,5);stats(congrus2,6)]);
plot_appear=false;
if plot_appear
    congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2);
    congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2);
    fh=figure('Color','w','Position',[100,100,100,235]);
    hold on
    swarmchart(ones(nnz(congrus1+congrus2),1),[stats(congrus1,5);stats(congrus2,6)],4,[0.8,0.8,0.8],'o','MarkerFaceAlpha',0.2,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','none')
    boxplot([stats(congrus1,5);stats(congrus2,6)],'Colors','k','Symbol','','Widths',0.8)
    set(gca(),'XTick',[],'YTick',0:0.2:1,'YTickLabel',0:20:100)
    xlabel('Real-shuffle')
    ylabel('Appeared trials (%)')
    exportgraphics(fh,'FC_appearance.pdf')
end
out_same=plot_one('same',metas,stats,idmap,ratiomap,false);
out_h2l=plot_one('h2l',metas,stats,idmap,ratiomap,false);
out_l2h=plot_one('l2h',metas,stats,idmap,ratiomap,false);

figure('Color','w','Position',[100,100,235,235]);
hold on;
for pi=0:3:6
    if pi==0,mm=out_same.mm;cic=out_same.cic;cie=out_same.cie;...
    elseif pi==3,mm=out_l2h.mm;cic=out_l2h.cic;cie=out_l2h.cie;...
    else,mm=out_h2l.mm;cic=out_h2l.cic;cie=out_h2l.cie;
    end
    ch=bar(pi+1,mm(1),'FaceColor','w','EdgeColor','k');
    eh=bar(pi+2,mm(2),'FaceColor','k','EdgeColor','k');
    errorbar(pi+(1:2),mm,[cic(1),cie(1)]-mm,[cic(2),cie(2)]-mm,'k.');
end
ylabel('F.C. selectivity  index')
text(max(xlim()),max(ylim()),sprintf('p=%.3f',ranksum([dists1;dists2],[dists1e;dists2e])),'HorizontalAlignment','right','VerticalAlignment','top');
ylim([-0.1,0.7]);
set(gca,'XTick',1.5:3:7.5,'XTickLabel',{'Within reg.','Low.->high.','High.->Low.'},'XTickLabelRotation',45)
legend([ch,eh],{'Correct trials','Error trials'});
keyboard()

end
function out=plot_one(type,metas,stats,idmap,ratiomap,plot_)

pvsst_ratio=arrayfun(@(x) ratiomap(char(idmap.ccfid2reg(x))),metas(:,6:7));
switch type
    case 'same'
        congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2) & metas(:,6)==metas(:,7);
        congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2) & metas(:,6)==metas(:,7);
        ftitle='Within region';
    case 'l2h'
        congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2) & pvsst_ratio(:,1)>pvsst_ratio(:,2);
        congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2) & pvsst_ratio(:,1)>pvsst_ratio(:,2);
        ftitle='Lower to higher';
    case 'h2l'
        congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2) & pvsst_ratio(:,1)<pvsst_ratio(:,2);
        congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2) & pvsst_ratio(:,1)<pvsst_ratio(:,2);
        ftitle='Higher to lower';
    case 'any'
        congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2);
        congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2);
        %     case 'forward'
        %         congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2);
        %         congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2);
        %     case 'backward'
        %         congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2);
        %         congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2);
end


dists1=arrayfun(@(x) norm(stats(x,1:2)-sum(stats(x,1:2),2)./2),find(congrus1));
negdir1=stats(congrus1,2)>stats(congrus1,1);
dists1(negdir1)=-dists1(negdir1);

dists2=arrayfun(@(x) norm(stats(x,1:2)-sum(stats(x,1:2),2)./2),find(congrus2));
negdir2=stats(congrus2,2)<stats(congrus2,1);
dists2(negdir2)=-dists2(negdir2);

dists1e=arrayfun(@(x) norm(stats(x,3:4)-sum(stats(x,3:4),2)./2),find(congrus1));
negdir1e=stats(congrus1,4)>stats(congrus1,3);
dists1e(negdir1e)=-dists1e(negdir1e);

dists2e=arrayfun(@(x) norm(stats(x,3:4)-sum(stats(x,3:4),2)./2),find(congrus2));
negdir2e=stats(congrus2,4)<stats(congrus2,3);
dists2e(negdir2e)=-dists2e(negdir2e);

mm=[mean([dists1;dists2]),mean([dists1e;dists2e])];
cic=bootci(1000,@(x) mean(x),[dists1;dists2]);
cie=bootci(1000,@(x) mean(x),[dists1e;dists2e]);
out.mm=mm;
out.cic=cic;
out.cie=cie;
if plot_
    
    
    fh=figure('Color','w','Position',[100,100,900,300]);
    subplot(1,3,1);
    hold on;
    ph1=scatter(stats(congrus1,1),stats(congrus1,2),9,'MarkerFaceColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
    ph2=scatter(stats(congrus2,1),stats(congrus2,2),9,'MarkerFaceColor','b','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
    plot([-10,50],[-10,50],'--k')
    xlim([-10,50])
    ylim([-10,50])
    title('Correct trials')
    legend([ph1,ph2],{'S1 congruent','S2 congruent'},'Location','northoutside');
    xlabel('S1 F.C. real-shuffle');
    ylabel('S2 F.C. real-shuffle');
    subplot(1,3,2);
    hold on;
    ph1=scatter(stats(congrus1,3),stats(congrus1,4),9,'MarkerFaceColor','m','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
    ph2=scatter(stats(congrus2,3),stats(congrus2,4),9,'MarkerFaceColor','c','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
    plot([-10,50],[-10,50],'--k')
    xlim([-10,50])
    ylim([-10,50])
    title('Error trials')
    legend([ph1,ph2],{'S1 congruent','S2 congruent'},'Location','northoutside');
    xlabel('S1 F.C. real-shuffle');
    ylabel('S2 F.C. real-shuffle');
    
    
    subplot(1,3,3);
    hold on;
    % ch=histogram([dists1;dists2],-6:0.5:6,'Normalization','probability');
    % eh=histogram([dists1e;dists2e],-6:0.5:6,'Normalization','probability');
    
    bar(1,mm(1),'FaceColor','w','EdgeColor','k')
    bar(2,mm(2),'FaceColor','k','EdgeColor','k')
    errorbar(1:2,mm,[cic(1),cie(1)]-mm,[cic(2),cie(2)]-mm,'k.');
    set(gca(),'XTick',1:2,'XTickLabel',{'Correct','Error'},'XTickLabelRotation',45)
    ylabel('F.C. selectivity  index')
    text(max(xlim()),max(ylim()),sprintf('p=%.3f',ranksum([dists1;dists2],[dists1e;dists2e])),'HorizontalAlignment','right','VerticalAlignment','top');
    ylim([-0.1,0.7]);
    sgtitle(ftitle);
end

end