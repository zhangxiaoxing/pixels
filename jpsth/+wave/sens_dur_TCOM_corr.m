function fh=sens_dur_TCOM_corr(fcom)
fh=figure('Color','w','Position',[16,180,1600,400]);
% TODDO: double check SU wave-type
t=tiledlayout(1,4);
t.Title.String='Olfaction/duration wave timing correlation';
[rs,rp]=deal([]);
for dd=["d3","d6"]
    nexttile
    hold on;

    inter_reg=intersect(...
        intersect(ephys.getGreyRegs('range','grey'),...
        fcom.(dd).collection(:,2)),...
        fcom.dur.collection(:,2));

    coord=nan(numel(inter_reg),3);
    regs=cell(numel(inter_reg),1);

    for ri=1:numel(inter_reg)
        sens_idx=(strcmp(fcom.(dd).collection(:,2),inter_reg(ri)));
        dur_idx=(strcmp(fcom.dur.collection(:,2),inter_reg(ri)));
        yy=fcom.dur.collection{dur_idx,1}./4;
        xx=fcom.(dd).collection{sens_idx,1}./4;
        coord(ri,:)=[xx,yy,1];
        regs(ri)=inter_reg(ri);
        scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(inter_reg{ri},'large_area',true),'MarkerEdgeColor','none');
        text(xx,yy,inter_reg{ri},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(inter_reg{ri},'large_area',true));
    end
%     coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    xlabel('Sensory wave COM')
    ylabel('Duration wave COM')
    [rs(end+1),ps]=corr(coord(:,1),coord(:,2),'type',"Spearman");
    [rp(end+1),pp]=corr(coord(:,1),coord(:,2),'type',"Pearson");

%     set(gca(),'XScale','log');
    text(max(xlim()),max(ylim()), ...
        {sprintf('r = %.3f, p = %.3f',rs,ps); ...
        sprintf('r = %.3f, p = %.3f',rp,pp)}, ...
        'HorizontalAlignment','right','VerticalAlignment','top');
    title(dd)

%     if opt.export
%         exportgraphics(fh,sprintf('per_region_TCOM_FRAC_%d.pdf',opt.delay));
%     end
end
nexttile()
bar(1:2,[rs(1),rs(2)],'FaceColor','k')
set(gca(),'YLim',[-1,1],'YTick',-1:0.5:1,'XTick',1:2,'XTickLabel',{'3s','6s'})
title('sens-dur-timing-spearman')

nexttile()
bar(1:2,[rp(1),rp(2)],'FaceColor','k')
set(gca(),'YLim',[-1,1],'YTick',-1:0.5:1,'XTick',1:2,'XTickLabel',{'3s','6s'})
title('sens-dur-timing-pearson')