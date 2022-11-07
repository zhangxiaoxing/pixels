function fh=sust_trans_bar_w_mix(sel_meta)
meta=ephys.util.load_meta('skip_stats',true);

% for fn=reshape(fieldnames(sel_meta),1,[])
%     meta.(fn{1})=sel_meta.(fn{1});
% end
all_olf_sel=ismember(sel_meta.wave_id,5:6);
sust_olf_sel=all(sel_meta.p_olf<0.05,2) & all_olf_sel;
trans_olf_sel=all_olf_sel & ~sust_olf_sel;

all_dur_sel=ismember(sel_meta.wave_id,7:8);
sust_dur_sel=all(sel_meta.p_dur<0.05,2) & all_dur_sel;
trans_dur_sel=all_dur_sel & ~sust_dur_sel;

all_mix_sel=ismember(sel_meta.wave_id,1:4);
sust_mix_sel=(all(sel_meta.p_mux<0.05,2)  ...
    |all(sel_meta.p_olf<0.05 | sel_meta.p_dur<0.05,2)) & all_mix_sel;
trans_mix_sel=all_mix_sel & ~sust_mix_sel;

tot=numel(sel_meta.wave_id);

[sust_olf_hat,sust_olf_ci]=binofit(nnz(sust_olf_sel),tot);
[trans_olf_hat,trans_olf_ci]=binofit(nnz(trans_olf_sel),tot);

[sust_dur_hat,sust_dur_ci]=binofit(nnz(sust_dur_sel),tot);
[trans_dur_hat,trans_dur_ci]=binofit(nnz(trans_dur_sel),tot);

[sust_mix_hat,sust_mix_ci]=binofit(nnz(sust_mix_sel),tot);
[trans_mix_hat,trans_mix_ci]=binofit(nnz(trans_mix_sel),tot);

fh=figure('Color','w');
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
% text(1,mean(ylim()),num2str(nnz(ismember(meta.mem_type,[1 3]))),'Rotation',90,'FontSize',10)
% text(2,mean(ylim()),num2str(nnz(ismember(meta.mem_type,[2 4]))),'Rotation',90,'FontSize',10)
% text(max(xlim()),max(ylim()),num2str(numel(meta.mem_type)),'HorizontalAlignment','right','VerticalAlignment','top')

end
