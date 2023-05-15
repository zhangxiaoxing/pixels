% No longer primary branch as of 2023.05.13. Data still checks out.

function plot_switch_cont_decoding(wrs_mux_meta,opt)
arguments
    wrs_mux_meta
    opt.type (1,:) char {mustBeMember(opt.type,{'olf','dur','mix'})}
end
switch opt.type
    case 'olf'
        dur_ids=5:6;
        olf_ids=5:6;
    case 'dur'
        dur_ids=7:8;
        olf_ids=7:8;
    case 'mix'
        dur_ids=1:4;
        olf_ids=1:4;
end

mix4durSWT=pct.pct_decoding_correct_error(wrs_mux_meta,dur_ids,'lblidx',8,'n_su',50,'switch_trial',true,'rpt',250);% dur
rptn=numel(mix4durSWT.dur.c_result_50su);
[~,~,pdur]=crosstab([zeros(rptn,1);ones(rptn,1);repmat(2,rptn,1)],...
    [mix4durSWT.dur.c_result_50su;mix4durSWT.dur.cs_result_50su;mix4durSWT.dur.e_result_50su]);

[swmm_dur,swci_dur]=binofit([nnz(mix4durSWT.dur.c_result_50su),nnz(mix4durSWT.dur.cs_result_50su),nnz(mix4durSWT.dur.e_result_50su)],...
    [numel(mix4durSWT.dur.c_result_50su),numel(mix4durSWT.dur.cs_result_50su),numel(mix4durSWT.dur.e_result_50su),]);

mix4olfSWT=pct.pct_decoding_correct_error(wrs_mux_meta,olf_ids,'lblidx',5,'n_su',50,'switch_trial',true,'rpt',250);% dur
[~,~,polf]=crosstab([zeros(rptn,1);ones(rptn,1);repmat(2,rptn,1)],...
    [mix4olfSWT.olf.c_result_50su;mix4olfSWT.olf.cs_result_50su;mix4olfSWT.olf.e_result_50su]);

[swmm_olf,swci_olf]=binofit([nnz(mix4olfSWT.olf.c_result_50su),nnz(mix4olfSWT.olf.cs_result_50su),nnz(mix4olfSWT.olf.e_result_50su)],...
    [numel(mix4olfSWT.olf.c_result_50su),numel(mix4olfSWT.olf.cs_result_50su),numel(mix4olfSWT.olf.e_result_50su),]);


figure()
hold on
bh=bar([swmm_olf;swmm_dur]);
errorbar([bh.XEndPoints],[bh.YEndPoints],...
    [swci_olf(1,1),swci_dur(1,1),swci_olf(2,1),swci_dur(2,1),swci_olf(3,1),swci_dur(3,1)]-[bh.YEndPoints],...
    [swci_olf(1,2),swci_dur(1,2),swci_olf(2,2),swci_dur(2,2),swci_olf(3,2),swci_dur(3,2)]-[bh.YEndPoints],'k.')
bh(1).FaceColor='w';
bh(2).FaceColor=[0.5,0.5,0.5];
bh(3).FaceColor='k';

ylim([0.45,1])
yline(0.5,'k--')

set(gca(),'XTick',1:2,'XTickLabel',{'Olf.','Dur.'},'YTick',0.5:0.25:1,'YTickLabel',50:25:100)
legend(bh,{'Correct cont.','Correct switch','Error trials'})
ylabel('Classification accuracy (%)')

title(sprintf('%s neuron, %.3f, %.3f',opt.type,polf,pdur))
end