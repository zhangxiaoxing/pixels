%% basic stats >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
global_init;
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
eff_meta=ephys.effect_size_meta();

sens_efsz=max(abs(eff_meta.cohen_d_olf),[],2);
sens_win=[min(sens_efsz)./2,prctile(sens_efsz,[20:20:100])];

dur_efsz=max(abs(eff_meta.cohen_d_dur),[],2);
dur_win=[min(dur_efsz)./2,prctile(dur_efsz,[20:20:100])];

pct_meta=pct.get_pct_meta(eff_meta,sens_win,dur_win);
com_map=wave.get_pct_com_map(pct_meta,'curve',true);

%ALT -> moved to new script
% wrs_mux_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',true,'merge_mux',true);
% wrs_mux_meta=ephys.get_wrs_mux_meta();
% com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true);


%% show case >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% mix same bin>>>>>>>>>>>>>>>>>>>>>>>>>
for bb=1:3
    [mxd,dd]=maxk(abs(eff_meta.cohen_d_dur(:,bb)),600);
    [mxo,oo]=maxk(abs(eff_meta.cohen_d_olf(:,bb)),600);

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
%>>>>mix, alternate bin>>>>>>>>>>>>>>>>>>>>>>>

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
% olf only ====================
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
% dur only ====================
reshape(find(any(out.p_mux(:,1:2)<0.002,2)),1,[]);

for ii=[3770,5818]%reshape(find(any(wrs_mux_meta.p_mux(:,1:2)<0.003,2) & median(wrs_mux_meta.class_fr,[2 3])>2.5 & median(wrs_mux_meta.class_fr,[2 3])<40),1,[])
    scfh=ephys.sens_dur_SC(ii,meta,'skip_raster',false);%
    if ~isempty(scfh)
        sgtitle(scfh,"Cross bin, SU #"+num2str(ii)+", dur");
%         exportgraphics(scfh,sprintf('SC/SCDUR%5d.png',ii),'ContentType','image');
        
    end
end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% wave heatmap >>>>>>>>>>>>>>>>>>
%prime grade
% fh=wave.plot_pct_wave(eff_meta,'sort_by',6);
% sgtitle(gcf(),'multi, sort by 6s')
% fh=wave.plot_pct_wave(eff_meta,'comb_set',2,'sort_by',6);
% sgtitle(gcf(),'single-mod, sort by 6s')
% % lesser grade
% fh=wave.plot_pct_wave(eff_meta,'sort_by',6,'lesser_grade',true);
% sgtitle(gcf(),'multi, sort by 6s, lesser grade')
% fh=wave.plot_pct_wave(eff_meta,'comb_set',2,'sort_by',6,'lesser_grade',true);
% sgtitle(gcf(),'single-mod, sort by 6s, lesser grade')
% lesser grade

%% TODO wave-half-half

fh=wave.plot_pct_wave(pct_meta,com_map,'sort_by',6,'scale',[0,1]);
sgtitle(gcf(),'multi, sort by 6s, extended')
fh=wave.plot_pct_wave(pct_meta,com_map,'comb_set',2,'sort_by',6,'scale',[0,1]);
sgtitle(gcf(),'single-mod, sort by 6s, extended')

%ALT

fh=wave.plot_pct_wave(wrs_mux_meta,com_map,'sort_by',6,'scale',[0,1]);
sgtitle(gcf(),'multi, sort by 6s, expanded')
fh=wave.plot_pct_wave(wrs_mux_meta,com_map,'comb_set',2,'sort_by',6,'scale',[0,1]);
sgtitle(gcf(),'single-mod, sort by 6s, expanded')



% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% effect size distribution histogram colorbar
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
[fh,olf_dec_olf]=wave.pct_decoding(sens_efsz,sens_win,'n_su',[10,50,100,200,300,500],'lblidx',5,'cmap','parula','new_data',true,'calc_dec',true,'rpt',100);
[fh,dur_dec_dur]=wave.pct_decoding(dur_efsz,dur_win,'n_su',[10,50,100,200,300,500],'lblidx',8,'cmap','cool','new_data',true,'calc_dec',true,'rpt',100);

%% cross decoding
wave.pct_decoding(sens_efsz,sens_win,'n_su',[10,50,100,200,300,500],'lblidx',8,'cmap','parula','cross',true,'new_data',true,'calc_dec',true)
title('Rank by odor-encoding, decoding duration')
exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);
wave.pct_decoding(dur_efsz,dur_win,'n_su',[10,50,100,200,300,500],'lblidx',5,'cmap','cool','cross',true,'new_data',true,'calc_dec',true)
title('Rank by duration-encoding, decoding odor')
exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);

%% correct error decoding >>>>>>>>>>>>>>>>>
if false
    odorMulti=pct.pct_decoding_correct_error(pct_meta,1:4,'lblidx',5,'n_su',50);% odor
    odorSingle=pct.pct_decoding_correct_error(pct_meta,5:6,'lblidx',5,'n_su',50);% odor
    durMulti=pct.pct_decoding_correct_error(pct_meta,1:4,'lblidx',8,'n_su',50);% odor
    durSingle=pct.pct_decoding_correct_error(pct_meta,7:8,'lblidx',8,'n_su',50);% odor
    save('corr_err_pct_decoding.mat','odorMulti','odorSingle','durMulti','durSingle');
else
    load('corr_err_pct_decoding.mat','odorMulti','odorSingle','durMulti','durSingle');
end

% plot correct error decoding, maybe encapsulate >>>>>>>>>>>>>>>>>>>>>
figure('Color','w','Position',[100,100,275,235]);
tiledlayout(1,2)
nexttile();
hold on
mm=[mean(odorMulti.olf.c_result_50su),mean(odorMulti.olf.e_result_50su),...
    mean(odorSingle.olf.c_result_50su),mean(odorSingle.olf.e_result_50su)];
bh=bar([1:2,4:5],diag(mm),'stacked');
[bh(1).FaceColor,bh(3).FaceColor]=deal('w');
[bh(2).FaceColor,bh(4).FaceColor]=deal('k');
sem=[std(odorMulti.olf.c_result_50su),std(odorMulti.olf.e_result_50su),...
    std(odorSingle.olf.c_result_50su),std(odorSingle.olf.e_result_50su)]...
    ./sqrt(numel(odorMulti.olf.c_result_50su));
errorbar([1:2,4:5],mm,sem,'k.');
ylim([0.4,1]);
set(gca(),'XTick',[1.5,4.5],'XTickLabel',{'Multi','Single'},'YTickLabel',get(gca(),'YTick').*100)
yline(0.5,'k--')
ylabel('Classification accuracy');

nexttile();
hold on
mm=[mean(durMulti.dur.c_result_50su),mean(durMulti.dur.e_result_50su),...
    mean(durSingle.dur.c_result_50su),mean(durSingle.dur.e_result_50su)];
bh=bar([1:2,4:5],diag(mm),'stacked');
[bh(1).FaceColor,bh(3).FaceColor]=deal('w');
[bh(2).FaceColor,bh(4).FaceColor]=deal('k');
sem=[std(durMulti.dur.c_result_50su),std(durMulti.dur.e_result_50su),...
    std(durSingle.dur.c_result_50su),std(durSingle.dur.e_result_50su)]...
    ./sqrt(numel(odorMulti.olf.c_result_50su));
errorbar([1:2,4:5],mm,sem,'k.');
ylim([0.4,1])
set(gca(),'XTick',[1.5,4.5],'XTickLabel',{'Multi','Single'},'YTickLabel',get(gca(),'YTick').*100)
yline(0.5,'k--')
ylabel('Classification accuracy');

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% mixed coding proportion matrix heatmap
win_cnt=numel(sens_win)-1;
count_mat=nan(win_cnt,win_cnt);
for s_idx=1:win_cnt
    for d_idx=1:win_cnt
        count_mat(s_idx,d_idx)=nnz(...
            pct_meta.mat_coord(:,1)==s_idx ...
            & pct_meta.mat_coord(:,2)==d_idx);
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
xlabel('Duration coding rank, higher is better (%)');
ylabel('Olfactory coding rank, higher is better (%)');
truesize(gcf(),[640,640])
% exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);


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

[map_cells,pct_bar_fh]=ephys.pct_reg_bars(pct_meta); % only need map_cells for tcom-frac corr
mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(map_cells,corr_log_log, ...
    'range','grey','data_type','pct-frac','stats_type','percentile',...
    'feat_tag',{'Mixed','Olfactory','Duration'},'corr2',true,'plot2',true);


%% TCOM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pct_tcom_fh=struct();
tcom_maps=cell(1,3);

% map_cells: mixed_map,olf_map,dur_map
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

mix_TCOM_GLM_2F_fh=wave.connectivity_proportion_GLM(tcom_maps(1),corr_ln_log, ...
    'range','grey','data_type','pct-TCOM','stats_type','percentile',...
    'feat_tag',{'Mixed'},'corr2',true,'plot2',true,'corr1',false);


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

%% TODO: COM_CHAIN
% K:\code\jpsth\+wave\COM_chain_SC.m

[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
bz.conn_prob_spatial_dist(sig,pair);
fc.fc_com_reg_wave
%>>> jump to TCOM section as needed
fh4=bz.inter_wave_pct(pct_meta);

%% TODO: FC_Decoding
% K:\code\jpsth\+bz\+fccoding\plot_coding.m

%% FC_TCOM_hierachy
% [fc_com_pvsst_stats_mix,fh_mix]=pct.fc_com_pct(com_map,pct_meta,'pair_type','congru');
% 
% [fc_com_pvsst_stats_mix,fh_mix]=pct.fc_com_pct(com_map,pct_meta,'pair_type','congru');
% 
% [fc_com_pvsst_stats_mix,fh_mix]=pct.fc_com_pct(com_map,pct_meta,'pair_type','incong');


bz.fc_conn_screen(com_map,pct_meta,'title_suffix','expanded')
%% exports

close all
% savefig(fhandles,sprintf('Ranksum1%s.fig',gather_config.fnsuffix));

