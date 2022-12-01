function plot_switch_cont_decoding(wrs_mux_meta,opt)
arguments
    wrs_mux_meta
    opt.merge_both (1,1) logical = true
end
if opt.merge_both
    dur_ids=[1:4,7:8];
    olf_ids=1:6;
else
    dur_ids=7:8;
    olf_ids=5:6;
end

mix4durSWT=pct.pct_decoding_correct_error(wrs_mux_meta,dur_ids,'lblidx',8,'n_su',50,'switch_trial',true,'rpt',250);% dur
rptn=numel(mix4durSWT.dur.c_result_50su);
[~,~,pdur]=crosstab((1:(2*rptn))>rptn,...
    [mix4durSWT.dur.c_result_50su;mix4durSWT.dur.cs_result_50su]);

[swmm_dur,swci_dur]=binofit([nnz(mix4durSWT.dur.c_result_50su),nnz(mix4durSWT.dur.cs_result_50su)],...
    [numel(mix4durSWT.dur.cs_result_50su),numel(mix4durSWT.dur.c_result_50su)]);

mix2olfSWT=pct.pct_decoding_correct_error(wrs_mux_meta,olf_ids,'lblidx',5,'n_su',50,'switch_trial',true,'rpt',250);% dur
[~,~,polf]=crosstab((1:(2*rptn))>rptn,...
    [mix2olfSWT.olf.c_result_50su;mix2olfSWT.olf.cs_result_50su]);

[swmm_olf,swci_olf]=binofit([nnz(mix2olfSWT.olf.c_result_50su),nnz(mix2olfSWT.olf.cs_result_50su)],...
    [numel(mix2olfSWT.olf.cs_result_50su),numel(mix2olfSWT.olf.c_result_50su)]);


figure()
hold on
bh=bar([swmm_olf;swmm_dur]);
errorbar([bh.XEndPoints],[bh.YEndPoints],...
    [swci_olf(1,1),swci_dur(1,1),swci_olf(2,1),swci_dur(2,1)]-[bh.YEndPoints],...
    [swci_olf(1,2),swci_dur(1,2),swci_olf(2,2),swci_dur(2,2)]-[bh.YEndPoints],'k.')
bh(1).FaceColor='w';
bh(2).FaceColor='k';

ylim([0.5,1])

set(gca(),'XTick',1:2,'XTickLabel',{'Olf.','Dur.'},'YTick',0.5:0.25:1,'YTickLabel',50:25:100)
legend(bh,{'Continuing trials','Switch trials'})
ylabel('Classification accuracy')
if opt.merge_both
    title(sprintf('with both %.3f, %.3f',polf,pdur))
else
    title(sprintf('without both %.3f, %.3f',polf,pdur))
end
end