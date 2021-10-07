function [fh,stats]=plot_stp_curve_correct_error(correct_stats,error_stats,opt)
%bz.hist.plot_stp_curve_correct_error(correct_stats,struct())
%[correct_stats,correct_type]=bz.hist.get_stp_stats(6000,'BZWT','trialtype','correct')
%[error_stats,error_type]=bz.hist.get_stp_stats(6000,'BZWT','trialtype','error')
arguments
    correct_stats (1,1) struct
    error_stats (1,1) struct
    opt.levels (1,:) double = 5
    opt.plot_dual (1,1) logical = false
    opt.hier_anova (1,1) logical = false
    opt.plot_progres_regres (1,1) logical = false
    opt.plot_loop_noloop (1,1) logical = false
end
fh=figure('Color','w','Position',[100,100,750,255]);
fidx=1;
stats=struct();
for level=opt.levels
    if opt.plot_loop_noloop
        figure('Color','w','Position',[100,100,1200,400]);
        midx=1;
        for mtype=["congru","incong","nonmem"]
            subplot(1,3,midx);midx=midx+1;
            hold on;
            crloop=bz.hist.get_stats_by_mem_type(correct_stats,'same_loop',level);
            crnoloop=bz.hist.get_stats_by_mem_type(correct_stats,'no_loop',level);
            ph=plot(mean(fliplr(crloop.(mtype)(:,2:end))).*100,'r-');
            rh=plot(mean(fliplr(crnoloop.(mtype)(:,2:end))).*100,'b-');
            set(gca,'XTick',1.5:1:10.5,'XTickLabel',200:200:2000);
            legend([ph,rh],{'Same loop','No loop'});
            xlim([-0.5,10.5])
            xlabel('Time lag (ms)')
            ylabel('Post spike gain (%)')
            title(mtype);
        end
        keyboard()
    end
    
    
    if opt.plot_progres_regres
        figure('Color','w','Position',[100,100,800,400]);
        for prog_diff=[true,false]
            if prog_diff
                subplot(1,2,2);
                title('Cross region')
            else
                subplot(1,2,1);
                title('All F.C.')
            end
            hold on;
            %diff
            crprog=bz.hist.get_stats_by_mem_type(correct_stats,'progres',level,'nonmem',false,'incong',false,'prog_cross_reg',prog_diff);
            crregr=bz.hist.get_stats_by_mem_type(correct_stats,'regres',level,'nonmem',false,'incong',false,'prog_cross_reg',prog_diff);
            ph=plot(mean(fliplr(crprog.congru(:,2:end))).*100,'r-');
            rh=plot(mean(fliplr(crregr.congru(:,2:end))).*100,'b-');
            set(gca,'XTick',1.5:1:10.5,'XTickLabel',200:200:2000);
            legend([ph,rh],{'Progressive','Regressive'});
            xlim([-0.5,10.5])
            xlabel('Time lag (ms)')
            ylabel('Post spike gain (%)')
        end
        keyboard()
    end
    
    
    
    %same
    cr=bz.hist.get_stats_by_mem_type(correct_stats,'within',level);
    stats.same=cr;
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'within',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),2,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(cr,er,'Within regions',...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
    
    y=[reshape(cr.congru(:,2:end).',[],1);reshape(cr.incong(:,2:end).',[],1);reshape(cr.nonmem(:,2:end).',[],1)];
    dirg=[ones(size(cr.congru,1).*10,1);2*ones(size(cr.incong,1).*10,1);3*ones(size(cr.nonmem,1).*10,1)];
    bing=repmat((1:10).',size(cr.congru,1)+size(cr.incong,1)+size(cr.nonmem,1),1);
    anovan(y,{dirg,bing})
    
    
    %diff
    crdiff=bz.hist.get_stats_by_mem_type(correct_stats,'between',level);
    stats.same=crdiff;
    if opt.plot_dual
        er=bz.hist.get_stats_by_mem_type(error_stats,'between',level);
    else
        er=struct();
    end
    subplot(numel(opt.levels),2,fidx);fidx=fidx+1;
    bz.hist.plot_one_stp_curve(crdiff,er,'cross regions',...
        'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
    
    y=[reshape(crdiff.congru(:,2:end).',[],1);reshape(crdiff.incong(:,2:end).',[],1);reshape(crdiff.nonmem(:,2:end).',[],1)];
    dirg=[ones(size(crdiff.congru,1).*10,1);2*ones(size(crdiff.incong,1).*10,1);3*ones(size(crdiff.nonmem,1).*10,1)];
    bing=repmat((1:10).',size(crdiff.congru,1)+size(crdiff.incong,1)+size(crdiff.nonmem,1),1);
    anovan(y,{dirg,bing})
    keyboard()
%     exportgraphics(fh,'STF_hier.pdf');
    
    %L2H vs H2L
    crl2h=bz.hist.get_stats_by_mem_type(correct_stats,'Low2High',level);
    crh2l=bz.hist.get_stats_by_mem_type(correct_stats,'High2Low',level);
    
    l2hmm=mean(fliplr(crl2h.congru(:,2:end)));
    l2hci=bootci(500,@(x) mean(x),fliplr(crl2h.congru(:,2:end)));
    
    h2lmm=mean(fliplr(crh2l.congru(:,2:end)));
    h2lci=bootci(500,@(x) mean(x),fliplr(crh2l.congru(:,2:end)));
    
    fh=figure('Color','w','Position',[32,32,200,200]);
    hold on;
    fill([200:200:2000,fliplr(200:200:2000)],[l2hci(1,:),fliplr(l2hci(2,:))],'r','EdgeColor','none','FaceAlpha',0.2);
    fill([200:200:2000,fliplr(200:200:2000)],[h2lci(1,:),fliplr(h2lci(2,:))],'m','EdgeColor','none','FaceAlpha',0.2);
    l2hh=plot(200:200:2000,l2hmm,'-r');
    h2lh=plot(200:200:2000,h2lmm,'-m');
    ylim([-0,0.05]);
    xlim([0,1000])
    set(gca,'XTick',500:500:2000,'YTick',0:0.02:0.05,'YTickLabel',0:2:5)
    legend([l2hh,h2lh],{'S to M','M to S'},'Location','northoutside');
    xlabel('Time lag (ms)')
    ylabel('Post spike gain(%)')
    exportgraphics(fh,'STF_h2l_l2h.pdf')
    
    y=[reshape(crl2h.congru(:,7:end).',[],1);reshape(crh2l.congru(:,7:end).',[],1)];
    dirg=[ones(size(crl2h.congru,1).*5,1);2*ones(size(crh2l.congru,1).*5,1)];
    bing=repmat((1:5).',size(crl2h.congru,1)+size(crh2l.congru,1),1);
    anovan(y,{dirg,bing})
    
    pperbin=nan(1,10);
    for pi=2:11
        pperbin(pi-1)=ranksum(crh2l.congru(:,pi),crl2h.congru(:,pi));
    end
    
    
    
    %     %L2H
    %     crl2h=bz.hist.get_stats_by_mem_type(correct_stats,'Low2High',level);
    %     stats.l2h=crl2h;
    %     if opt.plot_dual
    %         er=bz.hist.get_stats_by_mem_type(error_stats,'Low2High',level);
    %     else
    %         er=struct();
    %     end
    %     subplot(numel(opt.levels),3,fidx);fidx=fidx+1;
    %     bz.hist.plot_one_stp_curve(crl2h,er,'Low to high',...
    %         'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
    %     %H2L
    %     crh2l=bz.hist.get_stats_by_mem_type(correct_stats,'High2Low',level);
    %     stats.h2l=crh2l;
    %     if opt.plot_dual
    %         er=bz.hist.get_stats_by_mem_type(error_stats,'High2Low',level);
    %     else
    %         er=struct();
    %     end
    %     subplot(numel(opt.levels),3,fidx);fidx=fidx+1;
    %     bz.hist.plot_one_stp_curve(crh2l,er,'High to low',...
    %         'plot_dual',opt.plot_dual,'cmp_label','error','ref_label','correct');
end

data={cr,crdiff};
for didx=1:2
    y=[reshape(data{didx}.congru(:,2:end).',[],1);reshape(data{didx}.nonmem(:,2:end).',[],1);reshape(data{didx}.incong(:,2:end).',[],1)];
    dirg=[zeros(size(data{didx}.congru,1).*10,1);ones(size(data{didx}.nonmem,1).*10,1);2*ones(size(data{didx}.incong,1).*10,1)];
    bing=repmat((1:10).',size(data{didx}.congru,1)+size(data{didx}.nonmem,1)+size(data{didx}.incong,1),1);
    anovan(y,{dirg,bing})
end
keyboard()
exportgraphics(fh,'STF_hier.pdf');
end
