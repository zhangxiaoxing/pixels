%% basic stats >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
keyboard()
global_init;
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);

% wrs_mux_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',true,'merge_mux',true,'extend6s',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false);


% map_cells: mixed_map,olf_map,dur_map
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
tcom3_maps=cell(1,3);
for typeidx=1:3
    type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
    [fcom3.(type).collection,fcom3.(type).com_meta]=wave.per_region_COM(...
        com_map,'sel_type',type,'com_field','com3');
    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom3.(type).collection(:,2));
    [~,tcidx]=ismember(ureg,fcom3.(type).collection(:,2));
    tcom3_maps{typeidx}=containers.Map(...
        ureg,num2cell(cellfun(@(x) x/4, fcom3.(type).collection(tcidx,1))));
end


tcom6_maps=cell(1,3);
for typeidx=1:3
    type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
    [fcom6.(type).collection,fcom6.(type).com_meta]=wave.per_region_COM(...
        com_map,'sel_type',type,'com_field','com6');
    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom6.(type).collection(:,2));
    [~,tcidx]=ismember(ureg,fcom6.(type).collection(:,2));
    tcom6_maps{typeidx}=containers.Map(...
        ureg,num2cell(cellfun(@(x) x/4, fcom6.(type).collection(tcidx,1))));
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

%% duration switch trial v continuation trial
ephys.plot_switch_cont_decoding(wrs_mux_meta,"type",'olf')
ephys.plot_switch_cont_decoding(wrs_mux_meta,"type",'dur')
ephys.plot_switch_cont_decoding(wrs_mux_meta,"type",'mix')

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%% Figure 2
fh=ephys.sust_trans_bar_w_mix(wrs_mux_meta);
fh=ephys.sust_trans_correct_error(wrs_mux_meta);



%% wave-half-half
if false
    rpt=100;
    com_halfs=cell(rpt,2);
    for ii=1:rpt
        [com_map_h1,com_map_h2]=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'rnd_half',true);
        com_halfs(ii,:)={com_map_h1,com_map_h2};
    end
    blame=vcs.blame();
    save(sprintf('com_halfs_%d.mat',rpt),'com_halfs','blame')
end
if false
    com_map_err=wave.get_pct_com_map(wrs_mux_meta,'curve',false,'err',true);
    blame=vcs.blame();
    save('com_error.mat','com_map_err','blame');
end
wave_half_half_fh=wave.plot_wave_half_half(sens_meta);
stats_half_half_fh=wave.COM_half_half(sens_meta);


%% wave 
for hbound=0.7
    fh=wave.plot_pct_wave(wrs_mux_meta,com_map,'flex_sort',true,'scale',[0,hbound],'gauss2d',true);
    sgtitle(gcf(),"multi, 0:"+num2str(hbound));
    fh=wave.plot_pct_wave(wrs_mux_meta,com_map,'comb_set',2,'flex_sort',true,'scale',[0,hbound],'gauss2d',true);
    sgtitle(gcf(),"single-mod 0:"+num2str(hbound))
end

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Proportion, TCOM related

[map_cells,pct_bar_fh]=ephys.pct_reg_bars(wrs_mux_meta,'xyscale',{'log','log'}); % only need map_cells for tcom-frac corr
[map_cells,pct_bar_fh]=ephys.pct_reg_bars(wrs_mux_meta,'xyscale',{'linear','linear'}); % only need map_cells for tcom-frac corr
ch=gcf().Children.Children;
ch(3).YLim=[0,0.6];

mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(map_cells,gather_config.corr_log_log, ...
    'range','grey','data_type','pct-frac','stats_type','percentile',...
    'feat_tag',{'Mixed','Olfactory','Duration'},'corr2',false,'plot2',false);


%% TCOM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if false
% 3s delay trials and 6s delay trials early 3s
mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(tcom_maps,gather_config.corr_ln_log, ...
    'range','grey','data_type','wrs-mux-TCOM','stats_type','wrs-mux','feat_tag',{'Mixed','Olfactory','Duration'});
end
% 3s delay trials only
mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(tcom3_maps,gather_config.corr_ln_log, ...
    'range','grey','data_type','wrs-mux-TCOM3','stats_type','wrs-mux3','feat_tag',{'Mixed','Olfactory','Duration'});

% 6s delay trials only
mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(tcom6_maps,gather_config.corr_ln_log, ...
    'range','grey','data_type','wrs-mux-TCOM6','stats_type','wrs-mux6','feat_tag',{'Mixed','Olfactory','Duration'});



olf_TCOM_GLM_2F_fh=wave.connectivity_proportion_GLM(tcom3_maps(2),gather_config.corr_ln_log, ...
    'range','grey','data_type','pct-TCOM','stats_type','percentile',...
    'feat_tag',{'Olfaction'},'corr2',true,'plot2',true,'corr1',false);

dur_TCOM_GLM_2F_fh=wave.connectivity_proportion_GLM(tcom3_maps(3),gather_config.corr_ln_log, ...
    'range','grey','data_type','pct-TCOM','stats_type','percentile',...
    'feat_tag',{'Duration'},'corr2',true,'plot2',true,'corr1',false);

% mix_TCOM_GLM_2F_fh=wave.connectivity_proportion_GLM(tcom_maps(1),gather_config.corr_ln_log, ...
%     'range','grey','data_type','pct-TCOM','stats_type','percentile',...
%     'feat_tag',{'Mixed'},'corr2',true,'plot2',true,'corr1',false);
% 

%% TCOM and proportion correlation
pct_tcom_fh3=struct();
pct_tcom_fh6=struct();
for typeidx=1:3
    type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
    conn_reg=subsref(["AON","AON","PIR"],struct(type='()',subs={{typeidx}}));
    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom3.(type).collection(:,2));
    ffrac.collection=...
        [num2cell(cellfun(@(x) x(1),map_cells{typeidx}.values(ureg))),...
        ureg,...
        num2cell(ones(numel(ureg),1)*5),...
        num2cell(cellfun(@(x) x(3),map_cells{typeidx}.values(ureg)))];

    [pct_tcom_fh3,t3]=wave.per_region_COM_frac(...
        fcom3.(type),ffrac,...
        'hier_reg',conn_reg,...
        'corr',gather_config.corr_type,...
        'sel_type',type);
    sgtitle(pct_tcom_fh3,type+" 3sec")
    
    [pct_tcom_fh6,t6]=wave.per_region_COM_frac(...
        fcom6.(type),ffrac,...
        'hier_reg',conn_reg,...
        'corr',gather_config.corr_type,...
        'sel_type',type);
    sgtitle(pct_tcom_fh6,type+" 6sec")

end


%% Functional coupling

%% TODO: COM_CHAIN
% K:\code\jpsth\+wave\COM_chain_SC.m


%% FIG 3 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%% FIG 4 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
% [sig,pair]=bz.load_sig_sums_conn_file('pair',true);
[fcstats3,fh6]=fc.fc_com_reg_wave(wrs_mux_meta,com_map,tcom6_maps,'delay',6);
[fcstats3,fh3]=fc.fc_com_reg_wave(wrs_mux_meta,com_map,tcom6_maps,'delay',3);

fc.wave_stay_disappear(wrs_mux_meta)

if false
    inter_wave_fh=bz.inter_wave_ext_bars(wrs_mux_meta);  % dur, olf vs Isocortex, Straitum and Midbrain
    %skipped for current manuscript
end
% bz.inter_wave_ext_bars()

% separate wave timing cdf
% wave.mix_single_wave_timing

%>>> jump to TCOM section as needed
fh4=bz.inter_wave_pct(wrs_mux_meta); %congru vs incongru vs nonmem
fh4.fig.Children.Subtitle.String='Excitatory';
fh4i=bz.inter_wave_pct(wrs_mux_meta,'inhibit',true);
fh4i.fig.Children.Subtitle.String='Inhibitory';
bz.conn_prob_spatial_dist(sig,pair);
%% FC_Decoding
bz.fccoding.plot_coding(wrs_mux_meta,'dtype','olf')
bz.fccoding.plot_coding(wrs_mux_meta,'dtype','dur')

%% FC_TCOM_hierachy
% [fc_com_pvsst_stats_mix,fh_mix]=pct.fc_com_pct(com_map,pct_meta,'pair_type','congru');
% 
% [fc_com_pvsst_stats_mix,fh_mix]=pct.fc_com_pct(com_map,pct_meta,'pair_type','congru');
% 
% [fc_com_pvsst_stats_mix,fh_mix]=pct.fc_com_pct(com_map,pct_meta,'pair_type','incong');


bz.fc_conn_screen(com_map,pct_meta,'title_suffix','expanded')

%% loops
% sums_all
bz.rings.ring_wave_freq(wrs_mux_meta); 
load(fullfile('bzdata','sums_ring_stats_all.mat'));
% bz.rings.rings_reg_pie(sums_all)
bz.rings.loop_occurance_per_reg_su(sums_all,su_meta);
bz.rings.rings_wave_dynamic(sums_all)
bz.rings.rings_su_wave_tcom_corr(sums_all)

%TODO: assembly time constant olf, both, 3s 6s

%% chain
wave.chain_stats;
wave.chains_time_constant
wave.chains_loops_sc
wave.chain_tag(chains) %statistics
wave.chain_SC %plot
wave.chain_sust_tag(chains,'burstInterval',300)
wave.chain_sust_tag(chains,'burstInterval',600)
wave.chain_sust_tag(chains,'burstInterval',1500)

%% exports



