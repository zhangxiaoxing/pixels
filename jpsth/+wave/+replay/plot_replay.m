function [fhb,fhs]=plot_replay(stats,xlbl,opt)
arguments
    stats double
    xlbl char
    opt.title (1,:) char = 'chains'
end

fhb=figure();
boxplot(stats.','Colors','k','Symbol','c.')
ylim([0,1.5])
set(gca(),'XTick',1:size(stats,1),'XTickLabel',xlbl,'YScale','linear')
for jj=2:size(stats,1)
    pp=ranksum(stats(1,:),stats(jj,:));
    text(jj,1.5,sprintf('%.3f',pp),'VerticalAlignment','top','HorizontalAlignment','center')
end
title(opt.title)
ylabel('Motif spike frequency (Hz)')

fhs=figure();
hold on
for ii=1:size(stats,1)
    swarmchart(repmat(ii,size(stats,2),1),stats(ii,:),16,'.')
end
ylim([-0.1,1.5])
set(gca(),'XTick',1:size(stats,1),'XTickLabel',xlbl,'YScale','linear','TickLabelInterpreter','none')
title(opt.title)
ylabel('Motif spike frequency (Hz)')
end

