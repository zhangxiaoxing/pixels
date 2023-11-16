global_init;
%%%%%%%%%% fraction
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

%%%%%%%%%%%% fraction vs anatomy
com_map=wave.get_pct_com_map(ln_sel_meta,'odor_only',true,'criteria','Learning');

fh3=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0,0.7],'gauss2d',true,'delay',3,'xlim',3);
fh6=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0,0.7],'gauss2d',true,'delay',6);

% fraction model
mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(ln_map,gather_config.corr_log_log, ...
    'range','grey','data_type','pct-frac','stats_type','percentile','criteria','Learning');


%%%%%%%%%%% TCOM vs anatomy

tcom3_maps=struct();
tcom6_maps=struct();
grey_regs=ephys.getGreyRegs('range','grey','criteria','Learning');

[fcom3.odor_only.collection,fcom3.odor_only.com_meta]=wave.per_region_COM(...
    com_map,'sel_type','odor_only','com_field','com3','criteria','Learning');
ureg=intersect(grey_regs,fcom3.odor_only.collection(:,2));
[~,tcidx]=ismember(ureg,fcom3.odor_only.collection(:,2));
tcom3_maps.odor_only=containers.Map(...
    ureg,num2cell(cellfun(@(x) x/4, fcom3.odor_only.collection(tcidx,1))));

[fcom6.odor_only.collection,fcom6.odor_only.com_meta]=wave.per_region_COM(...
    com_map,'sel_type','odor_only','com_field','com6','criteria','Learning');
ureg=intersect(grey_regs,fcom6.odor_only.collection(:,2));
[~,tcidx]=ismember(ureg,fcom6.odor_only.collection(:,2));
tcom6_maps.odor_only=containers.Map(...
    ureg,num2cell(cellfun(@(x) x/4, fcom6.odor_only.collection(tcidx,1))));


% TCOM model 3s delay trials only
mixed_TCOM3_GLM_fh=wave.connectivity_proportion_GLM(tcom3_maps,gather_config.corr_ln_log, ...
    'range','grey','data_type','wrs-mux-TCOM3','stats_type','wrs-mux3','criteria','Learning');

% TCOM model 6s delay trials only
mixed_TCOM6_GLM_fh=wave.connectivity_proportion_GLM(tcom6_maps,gather_config.corr_ln_log, ...
    'range','grey','data_type','wrs-mux-TCOM6','stats_type','wrs-mux6','criteria','Learning');