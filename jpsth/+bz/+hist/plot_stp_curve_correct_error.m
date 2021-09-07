function [fh,stats]=plot_stp_curve_correct_error(correct_stats,error_stats,opt)
 %bz.hist.plot_stp_curve_correct_error(correct_stats,struct())
%[correct_stats,correct_type]=bz.hist.get_stp_stats(6000,'BZWT','trialtype','correct')
%[error_stats,error_type]=bz.hist.get_stp_stats(6000,'BZWT','trialtype','error')
arguments
    correct_stats (1,1) struct
    error_stats (1,1) struct
    opt.levels (1,:) double = 5
    opt.plot_dual (1,1) logical = false
end
fh=figure('Color','w','Position',[100,100,750,255]);
fidx=1;
stats=struct();
for level=opt.levels
    %same    
    cr=bz.hist.get_stats_by_mem_type(correct_stats,'within',level);
    stats.same=cr;
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'within',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),3,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(cr,er,'Within regions',...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
    %L2H
    cr=bz.hist.get_stats_by_mem_type(correct_stats,'Low2High',level);
    stats.l2h=cr;
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'Low2High',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),3,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(cr,er,'Lower to higer',...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
    %H2L
    cr=bz.hist.get_stats_by_mem_type(correct_stats,'High2Low',level);
    stats.h2l=cr;
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'High2Low',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),3,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(cr,er,'Higher to lower',...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
end
keyboard()
exportgraphics(fh,'STF_hier.pdf');
end
