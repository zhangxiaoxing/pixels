wt_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","WT","load_file",false,"skip_stats",true);
wt_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','WT','extend6s',true);
ln_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Learning","load_file",false,"skip_stats",true);
ln_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);

[wt_map,wt_fh]=ephys.pct_reg_bars(wt_su_meta,wt_sel_meta,'xyscale',{'linear','linear'},'only_odor',true,'criteria','WT'); % only need map_cells for tcom-frac corr
[ln_map,ln_fh]=ephys.pct_reg_bars(ln_su_meta,ln_sel_meta,'xyscale',{'linear','linear'},'only_odor',true,'criteria','Learning'); % only need map_cells for tcom-frac corr


[sfrac,sidx]=sort(subsref(cell2mat(wt_map.olf.values.'),substruct('()',{':',1})),'descend');
regs=subsref(wt_map.olf.keys,substruct('()',{sidx}));
frac_mat=[sfrac,nan(size(sfrac))];
for lnkey=reshape(ln_map.olf.keys,1,[])
    [in,pos]=ismember(lnkey,regs);
    if in
        frac_mat(pos,2)=subsref(ln_map.olf(lnkey{1}),substruct('()',{1}));
    end
end

figure()
bar(frac_mat,'grouped')
set(gca,'XTick',1:numel(regs),'XTickLabel',regs);
