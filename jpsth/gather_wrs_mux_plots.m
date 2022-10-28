%% basic stats >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
global_init;
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);

% wrs_mux_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',true,'merge_mux',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true);

tcom_maps=cell(1,3);

% map_cells: mixed_map,olf_map,dur_map
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for typeidx=1:3
    type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
    [fcom.(type).collection,fcom.(type).com_meta]=wave.per_region_COM(...
        com_map,'sel_type',type);
    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom.(type).collection(:,2));
    [~,tcidx]=ismember(ureg,fcom.(type).collection(:,2));
    tcom_maps{typeidx}=containers.Map(ureg,num2cell(cellfun(@(x) x/4, fcom.(type).collection(tcidx,1))));
end


%% show case >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% olf >>>>>>>>>>>>>>>>>>>>>>>>>
% #510
% idx=find(ismember(wrs_mux_meta.wave_id,5:6) & all(wrs_mux_meta.p_olf<1e-12,2));
idx=[510];
for ii=reshape(idx,1,[])
    scfh=ephys.sens_dur_SC(ii,meta,'skip_raster',false,'skip_fill',true);%
    if ~isempty(scfh)
        sgtitle(scfh, "SU #"+num2str(ii)+", OLF");
%         keyboard();
    end
end

%<<<<<<<<<<<<<<<<<<<<<<
%>>>>dur >>>>>>>>>>>>>>>>>>>>>>>
% idx=find(ismember(wrs_mux_meta.wave_id,7:8) & any(wrs_mux_meta.p_dur<1e-4,2));
idx=[2617];
for ii=reshape(idx,1,[])
    scfh=ephys.sens_dur_SC(ii,meta,'skip_raster',false,'skip_fill',true);%
    if ~isempty(scfh)
        sgtitle(scfh, "SU #"+num2str(ii)+", DUR");
        if false
            exportgraphics(scfh,fullfile("SC","dur"+num2str(ii)+".png"));
            close(scfh);
        else
%             keyboard()
        end
    end
end

% mux ====================
% idx=find(ismember(wrs_mux_meta.wave_id,1:4) & any(wrs_mux_meta.p_mux(:,1:2)<1e-2,2));
idx=[5818 23639];
for ii=reshape(idx,1,[])
    scfh=ephys.sens_dur_SC(ii,meta,'skip_raster',false,'skip_fill',true);%
    if ~isempty(scfh)
        sgtitle(scfh,"Cross bin, SU #"+num2str(ii)+", mux");
        if false
            exportgraphics(scfh,fullfile("SC","mux"+num2str(ii)+".png"));
            close(scfh);
        else
%             keyboard()
        end
    end
end


% %% svm decoding
% [fh,olf_dec_olf]=wave.pct_decoding(sens_efsz,sens_win,'n_su',[10,50,100,200,300,500],'lblidx',5,'cmap','parula','new_data',true,'calc_dec',true,'rpt',100);
% [fh,dur_dec_dur]=wave.pct_decoding(dur_efsz,dur_win,'n_su',[10,50,100,200,300,500],'lblidx',8,'cmap','cool','new_data',true,'calc_dec',true,'rpt',100);
% 
% %% cross decoding
% wave.pct_decoding(sens_efsz,sens_win,'n_su',[10,50,100,200,300,500],'lblidx',8,'cmap','parula','cross',true,'new_data',true,'calc_dec',true)
% title('Rank by odor-encoding, decoding duration')
% exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);
% wave.pct_decoding(dur_efsz,dur_win,'n_su',[10,50,100,200,300,500],'lblidx',5,'cmap','cool','cross',true,'new_data',true,'calc_dec',true)
% title('Rank by duration-encoding, decoding odor')
% exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);

%% correct error decoding >>>>>>>>>>>>>>>>>
if false
    odor4odor=pct.pct_decoding_correct_error(wrs_mux_meta,5:6,'lblidx',5,'n_su',50);% odor
    dur4odor=pct.pct_decoding_correct_error(wrs_mux_meta,7:8,'lblidx',5,'n_su',50);% odor
    mux4odor=pct.pct_decoding_correct_error(wrs_mux_meta,1:4,'lblidx',5,'n_su',50);% odor

    odor4dur=pct.pct_decoding_correct_error(wrs_mux_meta,5:6,'lblidx',8,'n_su',50);% dur
    dur4dur=pct.pct_decoding_correct_error(wrs_mux_meta,7:8,'lblidx',8,'n_su',50);% dur
    mux4dur=pct.pct_decoding_correct_error(wrs_mux_meta,1:4,'lblidx',8,'n_su',50);% dur
    save('corr_err_wrs_mux_decoding.mat','odor4odor','dur4odor','mux4odor','odor4dur','dur4dur','mux4dur');
else
    load('corr_err_wrs_mux_decoding.mat','odor4odor','dur4odor','mux4odor','odor4dur','dur4dur','mux4dur');
end
fh=ephys.plot_decode_correct_error(odor4odor,odor4dur,dur4odor,dur4dur,mux4odor,mux4dur);


%% TODO wave-half-half
rpt=2;
com_halfs=cell(rpt,2);
for ii=1:rpt
    [com_map_h1,com_map_h2]=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'rnd_half',true);
    com_halfs(ii,:)={com_map_h1,com_map_h2};
end


%%
fh=wave.plot_pct_wave(wrs_mux_meta,com_map,'sort_by',6,'scale',[0,1]);
sgtitle(gcf(),'multi, sort by 6s, expanded')
fh=wave.plot_pct_wave(wrs_mux_meta,com_map,'comb_set',2,'sort_by',6,'scale',[0,1]);
sgtitle(gcf(),'single-mod, sort by 6s, expanded')

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Proportion, TCOM related

[map_cells,pct_bar_fh]=ephys.pct_reg_bars(wrs_mux_meta,'xyscale',{'log','log'}); % only need map_cells for tcom-frac corr
[map_cells,pct_bar_fh]=ephys.pct_reg_bars(wrs_mux_meta,'xyscale',{'linear','linear'}); % only need map_cells for tcom-frac corr
ch=gcf().Children.Children;
ch(3).YLim=[0,0.6];

mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(map_cells,corr_log_log, ...
    'range','grey','data_type','pct-frac','stats_type','percentile',...
    'feat_tag',{'Mixed','Olfactory','Duration'},'corr2',true,'plot2',true);


%% TCOM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(tcom_maps,corr_ln_log, ...
    'range','grey','data_type','wrs-mux-TCOM','stats_type','wrs-mux','feat_tag',{'Mixed','Olfactory','Duration'});

olf_TCOM_GLM_2F_fh=wave.connectivity_proportion_GLM(tcom_maps(2),corr_ln_log, ...
    'range','grey','data_type','pct-TCOM','stats_type','percentile',...
    'feat_tag',{'Olfaction'},'corr2',true,'plot2',true,'corr1',false);

dur_TCOM_GLM_2F_fh=wave.connectivity_proportion_GLM(tcom_maps(3),corr_ln_log, ...
    'range','grey','data_type','pct-TCOM','stats_type','percentile',...
    'feat_tag',{'Duration'},'corr2',true,'plot2',true,'corr1',false);

% mix_TCOM_GLM_2F_fh=wave.connectivity_proportion_GLM(tcom_maps(1),corr_ln_log, ...
%     'range','grey','data_type','pct-TCOM','stats_type','percentile',...
%     'feat_tag',{'Mixed'},'corr2',true,'plot2',true,'corr1',false);
% 

%% TCOM and proportion correlation
pct_tcom_fh=struct();
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
fc.fc_com_reg_wave(wrs_mux_meta,com_map,tcom_maps)
%>>> jump to TCOM section as needed
fh4=bz.inter_wave_pct(wrs_mux_meta);
fh4.fig.Children.Subtitle.String='Excitatory';
fh4i=bz.inter_wave_pct(wrs_mux_meta,'inhibit',true);
fh4i.fig.Children.Subtitle.String='Inhibitory';
bz.conn_prob_spatial_dist(sig,pair);
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

