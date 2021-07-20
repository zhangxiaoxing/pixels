function fh=plot_stp_curve_correct_error(correct_stats,error_stats,opt)
%[correct_stats,correct_type]=bz.hist.get_stp_stats(6000,'BZWT','trialtype','correct')
%[error_stats,error_type]=bz.hist.get_stp_stats(6000,'BZWT','trialtype','error')
arguments
    correct_stats (1,1) struct
    error_stats (1,1) struct
    opt.levels (1,:) double = 5
    opt.plot_dual (1,1) logical = false
end
fh=figure('Color','w');
fidx=1;

for level=opt.levels
    %diff
    cr=bz.hist.get_stats_by_mem_type(correct_stats,'between',level);
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'between',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),2,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(cr,er,sprintf('Between regions %d',level),...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
    %same    
    cr=bz.hist.get_stats_by_mem_type(correct_stats,'within',level);
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'within',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),2,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(cr,er,sprintf('Within regions %d',level),...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
end
end
