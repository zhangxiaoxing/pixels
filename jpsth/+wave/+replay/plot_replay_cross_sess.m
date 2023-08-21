function fhb=plot_replay_cross_sess(cross_sess_mat,xlbl,opt)
arguments
    cross_sess_mat cell
    xlbl char
    opt.title (1,:) char = 'chains'
    opt.median_value (1,1) logical = false
end
mergemat=cell2mat(arrayfun(@(x) [cross_sess_mat{x},repmat(x,numel(cross_sess_mat{x}),1)],(1:numel(cross_sess_mat)).','UniformOutput',false));
mergemat=mergemat(isfinite(mergemat(:,1)),:);

fhb=figure('Position',[100,100,600,400]);
for ii=1:size(cross_sess_mat,2)
    boxplot(mergemat(:,1),mergemat(:,2),'Colors','k','Symbol','c.')
end
set(gca(),'XTick',1:size(cross_sess_mat,2),'XTickLabel',xlbl,'XTickLabelRotation',90,'YScale','log')

if opt.median_value
    for jj=1:size(cross_sess_mat,2)
        mm=median(cross_sess_mat{jj}(isfinite(cross_sess_mat{jj})));
        text(jj,1,sprintf('%.1f',mm),'VerticalAlignment','bottom','HorizontalAlignment','center')
    end
end

ylim([0.002,300]);
%
% fhs=figure();
% hold on
% for ii=1:size(stats_all,1)
%     swarmchart(repmat(ii,size(stats_all,2),1),stats_all(ii,:),16,'.')
% end
% ylim([-0.1,1.5])
% set(gca(),'XTick',1:size(stats_all,1),'XTickLabel',xlbl,'YScale','linear','TickLabelInterpreter','none')
title(opt.title)
ylabel('Motif spike frequency (Hz)')
end

