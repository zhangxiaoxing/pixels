function fhb=plot_replay_3panel(per_sess_mat,xlbl,opt)
arguments
    per_sess_mat double
    xlbl
    opt.title (1,:) char = 'chains'
    opt.median_value (1,1) logical = false
    opt.stats char {mustBeMember(opt.stats,{'mean','median'})} = 'median'
end

cmatv=reshape(per_sess_mat,[],1);
cmatg=reshape(repmat(1:size(per_sess_mat,2),size(per_sess_mat,1),1),[],1);
if strcmp(opt.stats,'median')
    mm=nanmedian(per_sess_mat);
    ci=bootci(1000,@(x) nanmedian(x),per_sess_mat);
else
    mm=nanmean(per_sess_mat);
    ci=nanstd(per_sess_mat)./sqrt(sum(isfinite(per_sess_mat)));
    ci=[ci;-ci]+mm;
end


fhb=figure('Position',[100,100,600,400]);
tiledlayout(1,3)
for pnl=1:2:5
    nexttile()
    hold on
    bar(mm(pnl:pnl+1).*eye(2),'stacked','EdgeColor','k')
    errorbar(1:2,mm(pnl:pnl+1),ci(1,pnl:pnl+1)-mm(pnl:pnl+1),ci(2,pnl:pnl+1)-mm(pnl:pnl+1),'k.');
    set(gca(),'XTick',1.5,'XTickLabel',xlbl{pnl})
    ylabel('Motif spike frequency (Hz)')
end

end
