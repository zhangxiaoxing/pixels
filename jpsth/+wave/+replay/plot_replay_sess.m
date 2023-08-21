function fhb=plot_replay_sess(per_sess_mat,xlbl,opt)
arguments
    per_sess_mat double
    xlbl char
    opt.title (1,:) char = 'chains'
    opt.ref_line (1,1) logical = false
    opt.ref_p_value (1,1) logical = true
    opt.median_value (1,1) logical = false
end
cmatv=reshape(per_sess_mat,[],1);
cmatg=reshape(repmat(1:size(per_sess_mat,2),size(per_sess_mat,1),1),[],1);


fhb=figure('Position',[100,100,600,400]);
boxplot(cmatv(isfinite(cmatv))+eps,cmatg(isfinite(cmatv)),'Colors','k','Symbol','c.')
if opt.ref_line
    yline(median(per_sess_mat(:,1)),'--r');
end
set(gca(),'XTick',1:size(per_sess_mat,2),'XTickLabel',xlbl,'XTickLabelRotation',90,'YScale','log')
if opt.ref_p_value
    for jj=2:size(per_sess_mat,2)
        finisel=isfinite(per_sess_mat(:,jj));
        pp=signrank(per_sess_mat(finisel,1),per_sess_mat(finisel,jj));
        text(jj,0.01,sprintf('%.3f',pp),'VerticalAlignment','bottom','HorizontalAlignment','center')
    end
end

if opt.median_value
    for jj=1:size(per_sess_mat,2)
        mm=median(cmatv(cmatg==jj & isfinite(cmatv)));
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
