function plot_ln_wt_stp_curve(learn_stats,wt_stats,opt)
% [learn_stats,~]=bz.hist.util.get_stp_stats(6000,'BZLN','suffix','Learning','type','neupix','criteria','Learning','any',true);
% [wt_stats,~]=bz.hist.util.get_stp_stats(6000,'BZWT','suffix','','type','neupix','criteria','WT','any',true);

arguments
    learn_stats (1,1) struct
    wt_stats (1,1) struct
    opt.levels (1,:) double = 3:5
end
figure('Color','w')
fidx=1;
for level=opt.levels
    %diff
    ln=bz.hist.learn.get_stats_by_mem_type(learn_stats,'between',level);
    wt=bz.hist.learn.get_stats_by_mem_type(wt_stats,'between',level);
    subplot(numel(opt.levels),2,fidx);fidx=fidx+1;
    bz.hist.learn.plotOne(wt,ln,sprintf('Between regions %d',level));
    %same    
    ln=bz.hist.learn.get_stats_by_mem_type(learn_stats,'within',level);
    wt=bz.hist.learn.get_stats_by_mem_type(wt_stats,'within',level);
    subplot(numel(opt.levels),2,fidx);fidx=fidx+1;
    bz.hist.learn.plotOne(wt,ln,sprintf('Within regions %d',level));
end
end
