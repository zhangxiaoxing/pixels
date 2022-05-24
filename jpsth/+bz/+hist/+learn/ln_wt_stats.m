[learn_stats,~]=bz.hist.util.get_stp_stats(6000,'BZLN','suffix','Learning','type','neupix','criteria','Learning','any',true);
[wt_stats,~]=bz.hist.util.get_stp_stats(6000,'BZWT','suffix','','type','neupix','criteria','WT','any',true);

%level3 ->{CTXpl,STR,TH,HY}
%level5 ->{PIR,AI,ORB,ILA,etc}
% stats_all=stats;
bz.hist.learn.plot_ln_wt_stp_curve(learn_stats,wt_stats)
bz.hist.learn.plot_ln_wt_per_reg(learn_stats,wt_stats,'levels',3:5)



