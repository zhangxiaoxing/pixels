wt_sel_meta=ephys.get_wrs_mux_meta('load_file',true,'save_file',false,'criteria','WT','extend6s',true);
ln_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);

[fh_wt,count_wt]=bz.inter_wave_pct(wt_sel_meta,'odor_only',true); %congru vs incongru vs nonmem bar lot
[fh_ln,count_ln]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'criteria','Learning'); %congru vs incongru vs nonmem bar lot


for regtype=["same_count","diff_count"]
    fh=figure();
    tiledlayout(1,3)
    for seltype=["nonmem","incong","congru"]
        nexttile()
        hold on
        bh=bar(diag([count_ln.(regtype).(seltype+"_hat"),count_wt.(regtype).(seltype+"_hat")]),'stacked');
        errorbar(1:2,...
            [count_ln.(regtype).(seltype+"_hat"),count_wt.(regtype).(seltype+"_hat")],...
            [count_ln.(regtype).(seltype+"_sem"),count_wt.(regtype).(seltype+"_sem")],'k.');
        ylabel('SC occurrence rate(%)')
        set(gca,'XTick',1:2,'XTickLabel',{'LN','WT'},'YTickLabel',get(gca,'YTick').*100)
        [~,~,pp]=crosstab([zeros(count_ln.(regtype).("pair_"+seltype),1);ones(count_wt.(regtype).("pair_"+seltype),1)],...
            [1:count_ln.(regtype).("pair_"+seltype)>count_ln.(regtype).("sig_"+seltype),1:count_wt.(regtype).("pair_"+seltype)>count_wt.(regtype).("sig_"+seltype)]);
        if pp<0.001
            text(1.5,mean(ylim),'***','HorizontalAlignment','center')
        elseif pp>0.05
            text(1.5,mean(ylim),'ns','HorizontalAlignment','center')
        else
            text(1.5,mean(ylim),sprintf('%.4f',pp,'HorizontalAlignment','center'));
        end
        title(seltype);
    end
    currtag=replace(regtype,"_count","-region");
    sgtitle(currtag);
    savefig(fh,fullfile("binary","ln_wt_SC_occur_rate-"+currtag+".fig"));
end