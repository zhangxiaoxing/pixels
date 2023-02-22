% The first 3s is always included due to necessasity for
% duration/mixed statistics measurement. Optional extends to include latter
% 3s of 6s duration

function fh=sust_trans_bar_w_mix(sel_meta,opt)
arguments
    sel_meta
    opt.exclude_switched (1,1) logical = true
    opt.merge_3_6 (1,1) logical = true
end
% meta=ephys.util.load_meta('skip_stats',true);
if opt.exclude_switched
    [olf_sw,dur_sw,mux_sw]=ephys.get_switched(sel_meta);
    sw_sel=olf_sw|dur_sw|mux_sw;
    sel_meta.wave_id(sw_sel)=-1;
else
    [olf_sw,dur_sw,mux_sw]=deal(false);
end
% for fn=reshape(fieldnames(sel_meta),1,[])
%     meta.(fn{1})=sel_meta.(fn{1});
% end
all_olf_sel=ismember(sel_meta.wave_id,5:6);
sust_olf_sel=all(sel_meta.p_olf<0.05,2) & all_olf_sel;
trans_olf_sel=all_olf_sel & ~sust_olf_sel;
% extend to 6s
sust_olf_sel6=all(sel_meta.p_olf<0.05,2) & all(sel_meta.p_olf6<0.05,2) & all_olf_sel;
trans_olf_sel6=all_olf_sel & ~sust_olf_sel6;


all_dur_sel=ismember(sel_meta.wave_id,7:8);
sust_dur_sel=all(sel_meta.p_dur<0.05,2) & all_dur_sel;
trans_dur_sel=all_dur_sel & ~sust_dur_sel;

all_mix_sel=ismember(sel_meta.wave_id,1:4);
sust_mix_sel=(all(sel_meta.p_mux<0.05,2)  ...
    |all(sel_meta.p_olf<0.05 | sel_meta.p_dur<0.05,2)) & all_mix_sel;
trans_mix_sel=all_mix_sel & ~sust_mix_sel;

sust_mix_sel6=(all(sel_meta.p_mux<0.05,2)  ...
    |all(sel_meta.p_olf<0.05 | sel_meta.p_dur<0.05,2))...
    & all(sel_meta.p_olf6<0.05,2) & all_mix_sel;
trans_mix_sel6=all_mix_sel & ~sust_mix_sel6;

tot=numel(sel_meta.wave_id);

[sust_olf_hat,sust_olf_ci]=binofit(nnz(sust_olf_sel),tot);
[trans_olf_hat,trans_olf_ci]=binofit(nnz(trans_olf_sel),tot);

[sust_dur_hat,sust_dur_ci]=binofit(nnz(sust_dur_sel),tot);
[trans_dur_hat,trans_dur_ci]=binofit(nnz(trans_dur_sel),tot);

[sust_mix_hat,sust_mix_ci]=binofit(nnz(sust_mix_sel),tot);
[trans_mix_hat,trans_mix_ci]=binofit(nnz(trans_mix_sel),tot);


% 6s extension
[sust_olf_hat6,sust_olf_ci6]=binofit(nnz(sust_olf_sel6),tot);
[trans_olf_hat6,trans_olf_ci6]=binofit(nnz(trans_olf_sel6),tot);

[sust_mix_hat6,sust_mix_ci6]=binofit(nnz(sust_mix_sel6),tot);
[trans_mix_hat6,trans_mix_ci6]=binofit(nnz(trans_mix_sel6),tot);

if opt.merge_3_6
    fh=figure('Color','w');
    tiledlayout(1,2);
    nexttile();
    hold on
    bh=bar([sust_olf_hat,trans_olf_hat;sust_olf_hat6,trans_olf_hat6],1,'grouped');
    bh(1).FaceColor='k';
    bh(2).FaceColor='w';

    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,[sust_olf_ci(1),sust_olf_ci6(1)]-bh(1).YEndPoints,...
        [sust_olf_ci(2),sust_olf_ci6(2)]-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,[trans_olf_ci(1),trans_olf_ci6(1)]-bh(2).YEndPoints,...
        [trans_olf_ci(2),trans_olf_ci6(2)]-bh(2).YEndPoints,'k.');
    ylabel('Fraction of all neurons (%)')
    set(gca(),'XTick',1:2,'XTickLabel',{'OLF3','OLF6'},'FontSize',10,'YTick',0:0.05:0.25,'YTickLabel',0:5:25)
    title('Olfactory')

    nexttile()
    hold on
    bh=bar([sust_mix_hat,trans_mix_hat;sust_mix_hat6,trans_mix_hat6],1,'grouped');
    bh(1).FaceColor='k';
    bh(2).FaceColor='w';

    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,[sust_mix_ci(1),sust_mix_ci6(1)]-bh(1).YEndPoints,...
        [sust_mix_ci(2),sust_mix_ci6(2)]-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,[trans_mix_ci(1),trans_mix_ci6(1)]-bh(2).YEndPoints,...
        [trans_mix_ci(2),trans_mix_ci6(2)]-bh(2).YEndPoints,'k.');
    ylabel('Fraction of all neurons (%)')
    set(gca(),'XTick',1:2,'XTickLabel',{'MIX3.','MIX6'},'FontSize',10,'YTick',0:0.05:0.25,'YTickLabel',0:5:25)
    title('Both-selective')
    ylim([0,0.1])

    disp({'sust','trans','switch','sust6','trans6'})
    disp({"olf",[nnz(sust_olf_sel),nnz(trans_olf_sel),nnz(olf_sw),nnz(sust_olf_sel6),nnz(trans_olf_sel6)]})
    disp({"dur",[nnz(sust_dur_sel),nnz(trans_dur_sel),nnz(dur_sw)]})
    disp({"mix",[nnz(sust_mix_sel),nnz(trans_mix_sel),nnz(mux_sw),nnz(sust_mix_sel6),nnz(trans_mix_sel6)]})
else
    fh=figure('Color','w');
    tiledlayout(1,2);
    nexttile()
    hold on
    bh=bar([sust_olf_hat,trans_olf_hat;sust_dur_hat,trans_dur_hat;sust_mix_hat,trans_mix_hat],1,'grouped');
    bh(1).FaceColor='k';
    bh(2).FaceColor='w';

    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,[sust_olf_ci(1),sust_dur_ci(1),sust_mix_ci(1)]-bh(1).YEndPoints,...
        [sust_olf_ci(2),sust_dur_ci(2),sust_mix_ci(2)]-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,[trans_olf_ci(1),trans_dur_ci(1),trans_mix_ci(1)]-bh(2).YEndPoints,...
        [trans_olf_ci(2),trans_dur_ci(2),trans_mix_ci(2)]-bh(2).YEndPoints,'k.');
    ylabel('Fraction of all neurons (%)')
    set(gca(),'XTick',1:3,'XTickLabel',{'Olf.','Dur.','Mixed'},'FontSize',10,'YTick',0:0.05:0.25,'YTickLabel',0:5:25)
    title('First 3s')

    nexttile()
    hold on
    bh=bar([sust_olf_hat6,trans_olf_hat6;sust_mix_hat6,trans_mix_hat6],1,'grouped');
    bh(1).FaceColor='k';
    bh(2).FaceColor='w';

    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,[sust_olf_ci6(1),sust_mix_ci6(1)]-bh(1).YEndPoints,...
        [sust_olf_ci6(2),sust_mix_ci6(2)]-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,[trans_olf_ci6(1),trans_mix_ci6(1)]-bh(2).YEndPoints,...
        [trans_olf_ci6(2),trans_mix_ci6(2)]-bh(2).YEndPoints,'k.');
    ylabel('Fraction of all neurons (%)')
    set(gca(),'XTick',1:2,'XTickLabel',{'Olf.','Mixed'},'FontSize',10,'YTick',0:0.05:0.25,'YTickLabel',0:5:25)
    title('Extends to 6s')

    disp({'sust','trans','switch','sust6','trans6'})
    disp({"olf",[nnz(sust_olf_sel),nnz(trans_olf_sel),nnz(olf_sw),nnz(sust_olf_sel6),nnz(trans_olf_sel6)]})
    disp({"dur",[nnz(sust_dur_sel),nnz(trans_dur_sel),nnz(dur_sw)]})
    disp({"mix",[nnz(sust_mix_sel),nnz(trans_mix_sel),nnz(mux_sw),nnz(sust_mix_sel6),nnz(trans_mix_sel6)]})
end
end
