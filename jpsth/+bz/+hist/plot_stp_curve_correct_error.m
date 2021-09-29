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
    crl2h=bz.hist.get_stats_by_mem_type(correct_stats,'Low2High',level);
    stats.l2h=crl2h;
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'Low2High',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),3,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(crl2h,er,'Low to high',...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
    %H2L
    crh2l=bz.hist.get_stats_by_mem_type(correct_stats,'High2Low',level);
    stats.h2l=crh2l;
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'High2Low',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),3,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(crh2l,er,'High to low',...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
end
y=[reshape(crl2h.congru(:,2:end).',[],1);reshape(crh2l.congru(:,2:end).',[],1)];
dirg=[ones(size(crl2h.congru,1).*10,1);2*ones(size(crh2l.congru,1).*10,1)];
bing=repmat((1:10).',size(crl2h.congru,1)+size(crh2l.congru,1),1);
anovan(y,{dirg,bing},'model','interaction')

pperbin=nan(1,10);
for pi=2:11
    pperbin(pi-1)=ranksum(crh2l.congru(:,pi),crl2h.congru(:,pi));
end

% y=[reshape(cr.congru(:,2:end).',[],1);reshape(crl2h.congru(:,2:end).',[],1);reshape(crh2l.congru(:,2:end).',[],1)];
% dirg=[zeros(size(cr.congru,1).*10,1);ones(size(crl2h.congru,1).*10,1);2*ones(size(crh2l.congru,1).*10,1)];
% bing=repmat((1:10).',size(cr.congru,1)+size(crl2h.congru,1)+size(crh2l.congru,1),1);
% anovan(y,{dirg,bing})


keyboard()
exportgraphics(fh,'STF_hier.pdf');
end
