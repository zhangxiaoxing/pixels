function plot_ln_wt_per_reg(learn_stats,wt_stats,opt)
arguments
    learn_stats (1,1) struct
    wt_stats (1,1) struct
    opt.levels (1,:) double = 3:5
end
for lvl=opt.levels
    regs=unique(cellfun(@(x) x{lvl}, learn_stats.any.reg(learn_stats.any.diff_reg(:,lvl)),'UniformOutput',false));
    lnreg=cellfun(@(x) x{lvl}, learn_stats.any.reg(learn_stats.any.diff_reg(:,lvl)),'UniformOutput',false);
    wtreg=cellfun(@(x) x{lvl}, wt_stats.any.reg(wt_stats.any.diff_reg(:,lvl)),'UniformOutput',false);
    lnfn=learn_stats.any.postspk(learn_stats.any.diff_reg(:,lvl),2:end);
    wtfn=wt_stats.any.postspk(wt_stats.any.diff_reg(:,lvl),2:end);
    pos=[];
    for ii=1:size(regs)
        onereg=regs{ii};
        lnsel=strcmp(lnreg,onereg);
        lnint=mean(sum(lnfn(lnsel,:),2));
        
        wtsel=strcmp(wtreg,onereg);
        wtint=mean(sum(wtfn(wtsel,:),2));
        
        pos(ii,1:2)=[lnint,wtint]*100;
    end
    figure('Color','w','Position',[100,100,400,300]);
    hold on
    for ii=1:size(regs)
        scatter(pos(ii,1),pos(ii,2),36,'k','MarkerFaceAlpha',0.2,'MarkerFaceColor','k','MarkerEdgeColor','none');
        text(pos(ii,1),pos(ii,2),regs{ii},'FontSize',12)
    end
    xlim([0,12])
    ylim([0,12])
    plot([0:12],[0:12],':r');
    title('Between-region integrated')
    xlabel('Learning (%)')
    ylabel('Welltrained (%)')
end