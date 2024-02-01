if ~exist('count_wt','var') || ~exist('count_ln','var') || ~exist('count_nv','var')
    wt_sel_meta=ephys.get_wrs_mux_meta('load_file',true,'save_file',false,'criteria','WT','extend6s',true);
    ln_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);
    nv_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Naive','extend6s',true);


    [fh_wt,count_wt]=bz.inter_wave_pct(wt_sel_meta,'odor_only',true); %congru vs incongru vs nonmem bar lot
    [fh_ln,count_ln]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'criteria','Learning'); %congru vs incongru vs nonmem bar lot
    [fh_nv,count_nv]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'criteria','Naive'); %congru vs incongru vs nonmem bar lot
end

for regtype=["same_count","diff_count"]
    fh=figure();
    tiledlayout(1,3)
    for seltype=["nonmem","incong","congru"]
        nexttile()
        hold on
        bh=bar(diag([count_nv.(regtype).(seltype+"_hat"),...
            count_ln.(regtype).(seltype+"_hat"),...
            count_wt.(regtype).(seltype+"_hat")]).*100,'stacked');
        errorbar(1:3,...
            100.*[count_nv.(regtype).(seltype+"_hat"),count_ln.(regtype).(seltype+"_hat"),count_wt.(regtype).(seltype+"_hat")],...
            100.*[count_nv.(regtype).(seltype+"_sem"),count_ln.(regtype).(seltype+"_sem"),count_wt.(regtype).(seltype+"_sem")],'k.');
        ylabel('SC occurrence rate(%)')
        set(gca,'XTick',1:3,'XTickLabel',{'NV','LN','WT'})
        [~,~,pp]=crosstab([...
            zeros(count_nv.(regtype).("pair_"+seltype),1);...
            ones(count_ln.(regtype).("pair_"+seltype),1);...
            2.*ones(count_wt.(regtype).("pair_"+seltype),1)],...
            [1:count_nv.(regtype).("pair_"+seltype)>count_nv.(regtype).("sig_"+seltype),......
            1:count_ln.(regtype).("pair_"+seltype)>count_ln.(regtype).("sig_"+seltype),...
            1:count_wt.(regtype).("pair_"+seltype)>count_wt.(regtype).("sig_"+seltype)]);
        if pp<0.001
            text(2,mean(ylim),'***','HorizontalAlignment','center')
        elseif pp>0.05
            text(2,mean(ylim),'ns','HorizontalAlignment','center')
        else
            text(2,mean(ylim),sprintf('%.4f',pp),'HorizontalAlignment','center');
        end
        title(seltype);
        if contains(regtype,"same")
            ylim([0,1.5]);
        else
            ylim([0,0.6])
        end
    end
    currtag=replace(regtype,"_count","-region");
    sgtitle(currtag);
    savefig(fh,fullfile("binary","ln_wt_SC_occur_rate-"+currtag+".fig"));
end