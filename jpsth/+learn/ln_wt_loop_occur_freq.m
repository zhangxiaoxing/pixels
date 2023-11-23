function ln_wt_loop_occur_freq(opt)
arguments
    opt.poolsize (1,1) double = 4
    opt.rpts (1,1) double = 100
end
[fh_ln,zscores_ln,data_ln,shuf_ln]=bz.rings.ring_wave_freq(repeats=opt.rpts,denovo=true,odor_only=true,poolsize=opt.poolsize,criteria='Learning');
pause(30); % wait for parpool cycle. It's a matlab bug.
[fh_wt,zscores_wt,data_wt,shuf_wt]=bz.rings.ring_wave_freq(repeats=opt.rpts,denovo=true,odor_only=true,poolsize=opt.poolsize,criteria='WT');

fh=figure();
tiledlayout(1,3)
for memtype=["nonmem","incongru","congru"]
    nexttile()
    hold on
    mm=[mean(zscores_ln.(memtype)(isfinite(zscores_ln.(memtype))));...
        mean(zscores_wt.(memtype)(isfinite(zscores_wt.(memtype))))];
    stdd=[std(zscores_ln.(memtype)(isfinite(zscores_ln.(memtype))));...
        std(zscores_wt.(memtype)(isfinite(zscores_wt.(memtype))))];
    pp=ranksum(zscores_ln.(memtype)(isfinite(zscores_ln.(memtype))),...
        zscores_wt.(memtype)(isfinite(zscores_wt.(memtype))));
    bh=bar(diag(mm),'stacked');
    sem=stdd./sqrt(sum(isfinite([zscores_ln.(memtype);zscores_wt.(memtype)]),2));
    errorbar(1:2,mm,sem,'k.');
    ylim([-25,500]);
    if pp<0.001
        text(1.5,400,'***','HorizontalAlignment','center')
    elseif frac_mat(pp,7)>0.05
        text(1.5,400,'ns','HorizontalAlignment','center')
    else
        text(1.5,400,sprintf('%.4f',frac_mat(pp,7)),'HorizontalAlignment','center')
    end
    title(memtype)
    ylabel('Z-score')
    set(gca,'XTick',1:2,'XTickLabel',{'LN','WT'})
end
sgtitle('Norm. loops occurr. rate');
savefig(fh,fullfile('binary','ln_wt_loop_occur_freq.fig'))
end