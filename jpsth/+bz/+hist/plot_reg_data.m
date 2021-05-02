ftick=6000;
[stats,memtypes]=bz.hist.util.get_stp_stats(ftick);
binsize=ftick./30;

% for treedep=5
%     bz.hist.plot_cross_reg(stats,'dir','to','treedepth',treedep,'relation','diff','title','Cross Region Input','binw',binsize);
%     bz.hist.plot_cross_reg(stats,'dir','from','treedepth',treedep,'relation','diff','title','Cross Region Output','binw',binsize);
%     bz.hist.plot_cross_reg(stats,'treedepth',treedep,'relation','same','title','Same region','binw',binsize);
% end

for treedep=5
    bz.hist.plot_cross_reg(stats,'dir','to','treedepth',treedep,'relation','diff','title','Cross Region Input','binw',binsize,'plottype','bar');
    bz.hist.plot_cross_reg(stats,'dir','from','treedepth',treedep,'relation','diff','title','Cross Region Output','binw',binsize,'plottype','bar');
    bz.hist.plot_cross_reg(stats,'treedepth',treedep,'relation','same','title','Same region','binw',binsize,'plottype','bar');
end