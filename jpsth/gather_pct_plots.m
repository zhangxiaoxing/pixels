%% show case >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
% mix >>>>>>>>>>>>>>>>>>>>>>>>>
for bb=1:3
    [mxd,dd]=maxk(abs(eff_meta.cohen_d_dur(:,bb)),500);
    [mxo,oo]=maxk(abs(eff_meta.cohen_d_olf(:,bb)),500);

    idx=intersect(dd,oo);
    for ii=reshape(idx,1,[])
        scfh=ephys.sens_dur_SC(ii,meta,'skip_raster',false);%
        if ~isempty(scfh)
            sgtitle(scfh,"Bin #"+num2str(bb)+", SU #"+num2str(ii)+", mixed");
            keyboard();
        end
    end
end
%<<<<<<<<<<<<<<<<<<<<<<
%>>>>mix, cross time>>>>>>>>>>>>>>>>>>>>>>>

[mxd,dd]=maxk(max(abs(eff_meta.cohen_d_dur),[],2),1000);
[mxo,oo]=maxk(max(abs(eff_meta.cohen_d_olf),[],2),1000);

idx=intersect(dd,oo);
for ii=reshape(idx,1,[])
    scfh=ephys.sens_dur_SC(ii,meta,'skip_raster',false);%
    if ~isempty(scfh)
        sgtitle(scfh,"Cross bin, SU #"+num2str(ii)+", mixed");
        keyboard();
    end
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

[mxd,dd]=mink(max(abs(eff_meta.cohen_d_dur),[],2),500);
[mxo,oo]=maxk(max(abs(eff_meta.cohen_d_olf),[],2),500);
idx=intersect(oo,dd);
for ii=reshape(idx,1,[])
    scfh=ephys.sens_dur_SC(ii,meta,'skip_raster',false);%
    if ~isempty(scfh)
        sgtitle(scfh,"Cross bin, SU #"+num2str(ii)+", olf");
        keyboard();
    end
end



[mxd,dd]=maxk(max(abs(eff_meta.cohen_d_dur),[],2),1000);
[mxo,oo]=mink(max(abs(eff_meta.cohen_d_olf),[],2),1000);
idx=intersect(oo,dd);
for ii=reshape(idx,1,[])
    scfh=ephys.sens_dur_SC(ii,meta,'skip_raster',false);%
    if ~isempty(scfh)
        sgtitle(scfh,"Cross bin, SU #"+num2str(ii)+", dur");
        keyboard();
    end
end




%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%% basic stats >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
eff_meta=ephys.effect_size_meta();
%% wave
fh=wave.plot_pct_wave(eff_meta,'sort_by',3);
sgtitle(gcf(),'multi, sort by 3s')
fh=wave.plot_pct_wave(eff_meta,'comb_set',2,'sort_by',3);
sgtitle(gcf(),'single-mod, sort by 3s')

fh=wave.plot_pct_wave(eff_meta,'sort_by',6);
sgtitle(gcf(),'multi, sort by 6s')
fh=wave.plot_pct_wave(eff_meta,'comb_set',2,'sort_by',6);
sgtitle(gcf(),'single-mod, sort by 6s')

fh=wave.plot_pct_wave(eff_meta,'sort_by',3,'xlim',3);
sgtitle(gcf(),'multi, sort by 3s')
fh=wave.plot_pct_wave(eff_meta,'comb_set',2,'sort_by',3,'xlim',3);
sgtitle(gcf(),'single-mod, sort by 3s')



sens_efsz=max(abs(eff_meta.cohen_d_olf),[],2);
sens_win=[min(sens_efsz)./2,prctile(sens_efsz,[20:20:100])];

dur_efsz=max(abs(eff_meta.cohen_d_dur),[],2);
dur_win=[min(dur_efsz)./2,prctile(dur_efsz,[20:20:100])];

% pct_meta=pct.get_pct_meta(eff_meta,sens_efsz,sens_win,dur_efsz,dur_win);

pct_meta4=pct.get_pct_meta(eff_meta,sens_efsz,sens_win,dur_efsz,dur_win,'single_mod_thresh',4);
pct_meta5=pct.get_pct_meta(eff_meta,sens_efsz,sens_win,dur_efsz,dur_win,'single_mod_thresh',5);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%% p-distribution histogram
cmap=flip(colormap(parula(5)));
bh=ephys.cmap_histogram(sens_efsz,sens_win,cmap);
set(gca(),'YScale','linear','XScale','linear')
set(gcf(),'Color','w')
set(gcf(),'Position',[100,100,520,100])
xlabel('p-value (permutation-test)')
ylabel('Number of neurons')
title('Olfactory guided 5 band binning')
% xlim([min(sens_efsz),max(sens_efsz)])
xlim([0,0.6])
exportgraphics(bh,'pct_decoding.pdf','ContentType','vector','Append',true);

cmap=colormap(cool(5));
bh=ephys.cmap_histogram(dur_efsz,dur_win,cmap);
set(gca(),'YScale','linear','XScale','linear')
set(gcf(),'Color','w')
set(gcf(),'Position',[100,100,520,100])
xlabel('p-value (permutation-test)')
ylabel('Number of neurons')
title('Duration guided 5 band binning')
xlim([0,0.6])
exportgraphics(bh,'pct_decoding.pdf','ContentType','vector','Append',true);

%% svm decoding
[fh,olf_dec_olf]=veve.pct_decoding(sens_efsz,sens_win,'n_su',[10,50,100,200,300,500],'lblidx',5,'cmap','parula','new_data',true,'calc_dec',true,'rpt',100);
[fh,dur_dec_dur]=wave.pct_decoding(dur_efsz,dur_win,'n_su',[10,50,100,200,300,500],'lblidx',8,'cmap','cool','new_data',true,'calc_dec',true,'rpt',100);

%% cross decoding
wave.pct_decoding(sens_efsz,sens_win,'n_su',[10,50,100,200,300,500],'lblidx',8,'cmap','parula','cross',true,'new_data',true,'calc_dec',true)
title('Rank by odor-encoding, decoding duration')
exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);
wave.pct_decoding(dur_efsz,dur_win,'n_su',[10,50,100,200,300,500],'lblidx',5,'cmap','cool','cross',true,'new_data',true,'calc_dec',true)
title('Rank by duration-encoding, decoding odor')
exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);



%% mixed coding proportion heatmap
win_cnt=numel(sens_win)-1;
count_mat=nan(win_cnt,win_cnt);
for s_idx=1:win_cnt
    for d_idx=1:win_cnt
        count_mat(s_idx,d_idx)=...
            nnz(sens_efsz>sens_win(s_idx) & ...
            sens_efsz<=sens_win(s_idx+1) & ...
            dur_efsz>dur_win(d_idx) & ...
            dur_efsz<=dur_win(d_idx+1));
    end
end

frac_mat=count_mat./sum(count_mat,'all');

figure('Color','w','Position',[100,100,640,640])
imagesc(frac_mat,[0.02,0.06])
colormap(flip(colormap('gray')))
cbh=colorbar();
set(gca(),'YDir','normal','XTick',0.5:1:10.5,'XTickLabel',0:20:100,...
    'YTick',0.5:1:10.5,'YTickLabel',0:20:100);
cbh.Label.String='Proportion of total population (%)';
cbh.TickLabels=cbh.Ticks*100;
xlabel('Olfactory coding rank, higher is better(%)');
ylabel('Duration coding rank, higher is better (%)');
truesize(gcf(),[640,640])
exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);


% %% region-dist
% if false
%     pct_map_cells=cell(1,4);
%     grey_regs=ephys.getGreyRegs('range','grey');
%     for mii=1:4
%         pct_map_cells{mii}=containers.Map('KeyType','char','ValueType','any');
%         for rr=grey_regs
%             reg_sel=strcmp(meta.reg_tree(5,:),rr).';
%             pct_map_cells{mii}(rr{1})=nnz(reg_sel & pct_meta.wave_id==mii)./nnz(reg_sel);
%         end
%     end
%     fh=ephys.pct_proportion_GLM(pct_map_cells,'PearsonLogLog',...
%         'range','grey','data_type',sprintf('class-%d bins',bw),'stats_type','pct',...
%         'feat_tag',{'Class#1','Class#2','Class#3','Class#4'});
%     for ff=reshape(fieldnames(fh),1,[])
%         exportgraphics(fh.(ff{1}),'pct_decoding.pdf','ContentType','vector','Append',true);
%     end
% end

%% Proportion, TCOM related
[map_cells,pct_bar_fh]=ephys.pct_reg_bars(pct_meta4); % only need map_cells for tcom-frac corr
mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(map_cells,corr_log_log, ...
    'range','grey','data_type','pct-frac','stats_type','percentile','feat_tag',{'Mixed','Olfactory','Duration'});


%% TCOM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pct_tcom_fh=struct();
tcom_maps=cell(1,3);

% map_cells: mixed_map,olf_map,dur_map
% com-map >>>>>>>>>>>>>>>>>>>>>>>>>>
global_init;
com_map=wave.get_pct_com_map(pct_meta4,'curve',true);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for typeidx=1:3
    type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
    [fcom.(type).collection,fcom.(type).com_meta]=wave.per_region_COM(...
        com_map,'pct_type',type);
    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom.(type).collection(:,2));
    [~,tcidx]=ismember(ureg,fcom.(type).collection(:,2));
    tcom_maps{typeidx}=containers.Map(ureg,num2cell(cellfun(@(x) x/4, fcom.(type).collection(tcidx,1))));
end

mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(tcom_maps,corr_ln_log, ...
    'range','grey','data_type','pct-TCOM','stats_type','percentile','feat_tag',{'Mixed','Olfactory','Duration'});

dur_TCOM_GLM_2F_fh=wave.connectivity_proportion_GLM(tcom_maps(3),corr_ln_log, ...
    'range','grey','data_type','pct-TCOM','stats_type','percentile',...
    'feat_tag',{'Duration'},'corr2',true,'plot2',true,'corr1',false);


%% TCOM and proportion correlation

for typeidx=1:3
    type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
    conn_reg=subsref(["EPd","AON","SS"],struct(type='()',subs={{typeidx}}));
    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom.(type).collection(:,2));
    ffrac.collection=...
        [num2cell(cellfun(@(x) x(1),map_cells{typeidx}.values(ureg))),...
        ureg,...
        num2cell(ones(numel(ureg),1)*5),...
        num2cell(cellfun(@(x) x(3),map_cells{typeidx}.values(ureg)))];

    [pct_tcom_fh,t]=wave.per_region_COM_frac(...
        fcom.(type),ffrac,...
        'hier_reg',conn_reg,...
        'corr',gather_config.corr_type,...
        'sel_type',type);
end


%% Functional coupling
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
pct_meta=pct.get_pct_meta();

%>>> jump to TCOM section as needed
if true
    % initiate global variable
    [fh,~,th]=bz.inter_wave_pct(pct_meta4);
    [fh,~,th]=bz.inter_wave_pct(pct_meta5);
    title(th,sprintf('class of %d bins',bw));
    exportgraphics(fh.mat,'pct_decoding.pdf','ContentType','vector','Append',true);
    exportgraphics(fh.bar,'pct_decoding.pdf','ContentType','vector','Append',true);
end



%% FC_TCOM_hierachy
[fc_com_pvsst_stats_mix,fh_mix]=pct.fc_com_pct(com_map,pct_meta,'pair_type','congru');
[fc_com_pvsst_stats_mix,fh_mix]=pct.fc_com_pct(com_map,pct_meta,'pair_type','incong');

% [fc_com_pvsst_stats_mix,fh_mix]=wave.fc_com_pvsst(com_map,struct(),pct_meta,'pct_stats',true,'hiermap','CP','descend',false,'mem_type','mixed');
% [fc_com_pvsst_stats_olf,fh_olf]=wave.fc_com_pvsst(com_map,struct(),pct_meta,'pct_stats',true,'hiermap','MOB','descend',true,'mem_type','olf');



%% exports
fhandles=get(groot(),'Children');
for hc=reshape(fhandles,1,[])
    exportgraphics(hc,'pct_decoding.pdf','ContentType','vector','Append',true);
end
% savefig(fhandles,sprintf('Ranksum1%s.fig',gather_config.fnsuffix));



