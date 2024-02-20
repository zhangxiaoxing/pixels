if ~exist('count_wt','var') || ~exist('count_ln','var') || ~exist('count_nv','var')
    wt_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','WT');
    ln_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','Learning');
    nv_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','Naive');

    [sc_fh.wt.bwide,sc_count.wt.bwide]=bz.inter_wave_pct(wt_sel_meta,'odor_only',true); %congru vs incongru vs nonmem bar lot
    [sc_fh.ln.bwide,sc_count.ln.bwide]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'criteria','Learning'); %congru vs incongru vs nonmem bar lot
    [sc_fh.nv.bwide,sc_count.nv.bwide]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'criteria','Naive'); %congru vs incongru vs nonmem bar lot

    [sc_fh.wt.MO,sc_count.wt.MO]=bz.inter_wave_pct(wt_sel_meta,'odor_only',true,'per_region',true,'region','MO'); %congru vs incongru vs nonmem bar lot
    [sc_fh.ln.MO,sc_count.ln.MO]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'per_region',true,'region','MO','criteria','Learning'); %congru vs incongru vs nonmem bar lot
    [sc_fh.nv.MO,sc_count.nv.MO]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'per_region',true,'region','MO','criteria','Naive'); %congru vs incongru vs nonmem bar lot

    [sc_fh.wt.OLF,sc_count.wt.OLF]=bz.inter_wave_pct(wt_sel_meta,'odor_only',true,'per_region',true,'region','OLF'); %congru vs incongru vs nonmem bar lot
    [sc_fh.ln.OLF,sc_count.ln.OLF]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'per_region',true,'region','OLF','criteria','Learning'); %congru vs incongru vs nonmem bar lot
    [sc_fh.nv.OLF,sc_count.nv.OLF]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'per_region',true,'region','OLF','criteria','Naive'); %congru vs incongru vs nonmem bar lot

    [sc_fh.wt.AI,sc_count.wt.AI]=bz.inter_wave_pct(wt_sel_meta,'odor_only',true,'per_region',true,'region','AI'); %congru vs incongru vs nonmem bar lot
    [sc_fh.ln.MO,sc_count.ln.AI]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'per_region',true,'region','AI','criteria','Learning'); %congru vs incongru vs nonmem bar lot
    [sc_fh.nv.AI,sc_count.nv.AI]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'per_region',true,'region','AI','criteria','Naive'); %congru vs incongru vs nonmem bar lot

    [sc_fh.wt.ORB,sc_count.wt.ORB]=bz.inter_wave_pct(wt_sel_meta,'odor_only',true,'per_region',true,'region','ORB'); %congru vs incongru vs nonmem bar lot
    [sc_fh.ln.ORB,sc_count.ln.ORB]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'per_region',true,'region','ORB','criteria','Learning'); %congru vs incongru vs nonmem bar lot
    [sc_fh.nv.ORB,sc_count.nv.ORB]=bz.inter_wave_pct(ln_sel_meta,'odor_only',true,'per_region',true,'region','ORB','criteria','Naive'); %congru vs incongru vs nonmem bar lot


end

for one_reg=["OLF","MO","AI","ORB"]
    for regtype=["same_count","diff_count"]
        fh=figure();
        tiledlayout(1,3)
        for seltype=["nonmem","incong","congru"]
            nexttile()
            hold on
            bh=bar(diag([sc_count.nv.(one_reg).(regtype).(seltype+"_hat"),...
                sc_count.ln.(one_reg).(regtype).(seltype+"_hat"),...
                sc_count.wt.(one_reg).(regtype).(seltype+"_hat")]).*100,'stacked');
            errorbar(1:3,...
                100.*[sc_count.nv.(one_reg).(regtype).(seltype+"_hat"),sc_count.ln.(one_reg).(regtype).(seltype+"_hat"),sc_count.wt.(one_reg).(regtype).(seltype+"_hat")],...
                100.*[sc_count.nv.(one_reg).(regtype).(seltype+"_sem"),sc_count.ln.(one_reg).(regtype).(seltype+"_sem"),sc_count.wt.(one_reg).(regtype).(seltype+"_sem")],'k.');
            ylabel('SC occurrence rate(%)')
            set(gca,'XTick',1:3,'XTickLabel',{'NV','LN','WT'})
            [~,~,pp]=crosstab([...
                zeros(sc_count.nv.(one_reg).(regtype).("pair_"+seltype),1);...
                ones(sc_count.ln.(one_reg).(regtype).("pair_"+seltype),1);...
                2.*ones(sc_count.wt.(one_reg).(regtype).("pair_"+seltype),1)],...
                [1:sc_count.nv.(one_reg).(regtype).("pair_"+seltype)>sc_count.nv.(one_reg).(regtype).("sig_"+seltype),......
                1:sc_count.ln.(one_reg).(regtype).("pair_"+seltype)>sc_count.ln.(one_reg).(regtype).("sig_"+seltype),...
                1:sc_count.wt.(one_reg).(regtype).("pair_"+seltype)>sc_count.wt.(one_reg).(regtype).("sig_"+seltype)]);
            if pp<0.001
                text(2,mean(ylim),'***','HorizontalAlignment','center')
            elseif pp>0.05
                text(2,mean(ylim),'ns','HorizontalAlignment','center')
            else
                text(2,mean(ylim),sprintf('%.4f',pp),'HorizontalAlignment','center');
            end
            title(seltype);
            if contains(regtype,"same")
                ylim([0,2]);
            else
                ylim([0,0.75])
            end
        end
        currtag=replace(regtype,"_count","-region");
        sgtitle(one_reg+": "+currtag);
        % savefig(fh,fullfile("binary","ln_wt_SC_occur_rate-"+currtag+".fig"));
    end
end