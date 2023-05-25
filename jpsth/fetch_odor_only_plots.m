% 37.72 w/o FDR, 17.53 w/FDR
keyboard()

%% basic stats >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
global_init;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_meta=ephys.get_wrs_meta('fdr',true);

com_map=wave.get_olf_com_map(wrs_meta,'curve',true);

grey_regs=ephys.getGreyRegs('range','grey');
[fcom.olf.collection,fcom.olf.com_meta]=wave.per_region_COM(...
    com_map,'sel_type','olf','com_field','com');
ureg=intersect(grey_regs,fcom.olf.collection(:,2));
[~,tcidx]=ismember(ureg,fcom.olf.collection(:,2));
reg_com_maps.olf=containers.Map(...
    ureg,num2cell(cellfun(@(x) x/4, fcom.olf.collection(tcidx,1))));

% on=nnz(ismember(wrs_meta.wave_id,5:6));
% alln=numel(wrs_meta.wave_id);

%% wave & stay
wave_n_stay=nnz(ismember(wrs_meta.wave_id,5:6) & wrs_meta.p_olf(:,3)<0.05 & all(wrs_meta.p_olf(:,4:6)<0.05,2));

%% show case >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% olf >>>>>>>>>>>>>>>>>>>>>>>>>
% #510
% idx=find(ismember(wrs_mux_meta.wave_id,5:6) & all(wrs_mux_meta.p_olf<1e-12,2));
idx=[510];
for ii=reshape(idx,1,[])
    scfh=ephys.sens_dur_SC(ii,su_meta,'skip_raster',false,'skip_fill',true);%
    if ~isempty(scfh)
        sgtitle(scfh, "SU #"+num2str(ii)+", OLF");
%         keyboard();
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
% svm on neuron firing rates
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
    [~,~,p]=crosstab(1:2000>1000,[odor4odor.olf.c_result_50su;odor4odor.olf.e_result_50su])
    [~,~,p]=crosstab(1:2000>1000,[odor4dur.dur.c_result_50su;odor4dur.dur.e_result_50su])
    [~,~,p]=crosstab(1:2000>1000,[mux4odor.olf.c_result_50su;mux4odor.olf.e_result_50su])
    [~,~,p]=crosstab(1:2000>1000,[mux4dur.dur.c_result_50su;mux4dur.dur.e_result_50su])
    [~,~,p]=crosstab(1:2000>1000,[dur4odor.olf.c_result_50su;dur4odor.olf.e_result_50su])
    [~,~,p]=crosstab(1:2000>1000,[dur4dur.dur.c_result_50su;dur4dur.dur.e_result_50su])
end
fh=ephys.plot_decode_correct_error(odor4odor,odor4dur,dur4odor,dur4dur,mux4odor,mux4dur);

%% duration switch trial v continuation trial
ephys.plot_switch_cont_decoding(wrs_mux_meta,"type",'olf')


%% sust trans
% TODO olf only
fh=ephys.sust_trans_bar_w_mix(wrs_meta);
fh=ephys.sust_trans_correct_error(wrs_meta);

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

while true
    wave_half_half_fh=wave.plot_wave_half_half(wrs_mux_meta,'minr',0.90);
    if ~isempty(wave_half_half_fh)
        break;
    end
end
stats_half_half_fh=wave.COM_half_half_wrs_mux();

%% wave 
for hbound=0.8
    fh=wave.plot_pct_wave(wrs_meta,com_map,'comb_set',3,'scale',[0,hbound],'gauss2d',true);
    sgtitle(gcf(),"single-mod 0:"+num2str(hbound))
end

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Proportion, TCOM related
[map_cells,pct_bar_fh]=ephys.pct_reg_bars(wrs_meta,'xyscale',{'linear','linear'},'only_odor',true); % only need map_cells for tcom-frac corr
% ch=gcf().Children.Children;
% ch(3).YLim=[0,0.6];
map_cells=rmfield(map_cells,{'mixed','dur'});

mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(map_cells,gather_config.corr_log_log, ...
    'range','grey','data_type','pct-frac','stats_type','percentile');

%% TCOM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% joint 3s-6s trials
mixed_TCOM6_GLM_fh=wave.connectivity_proportion_GLM(reg_com_maps,gather_config.corr_ln_log, ...
    'range','grey','data_type','wrs-TCOM','stats_type','wrs');



%% TCOM and proportion correlation % 4 panel scatters
% pct_tcom_fh3=struct();
% pct_tcom_fh6=struct();
% for typeidx=1:3
%     type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
%     conn_reg=subsref(["AON","AON","PIR"],struct(type='()',subs={{typeidx}}));
%     ureg=intersect(ephys.getGreyRegs('range','grey'),...
%         fcom3.(type).collection(:,2));
%     ffrac.collection=...
%         [num2cell(cellfun(@(x) x(1),map_cells.(type).values(ureg))),...
%         ureg,...
%         num2cell(ones(numel(ureg),1)*5),...
%         num2cell(cellfun(@(x) x(3),map_cells.(type).values(ureg)))];
% 
%     [pct_tcom_fh3,t3]=wave.per_region_COM_frac(...
%         fcom3.(type),ffrac,...
%         'hier_reg',conn_reg,...
%         'corr',gather_config.corr_type,...
%         'sel_type',type);
%     sgtitle(pct_tcom_fh3,type+" 3sec")
%     
%     [pct_tcom_fh6,t6]=wave.per_region_COM_frac(...
%         fcom.(type),ffrac,...
%         'hier_reg',conn_reg,...
%         'corr',gather_config.corr_type,...
%         'sel_type',type);
%     sgtitle(pct_tcom_fh6,type+" 6sec")
% 
% end


%% Functional coupling
fh4=bz.inter_wave_pct(wrs_meta); %congru vs incongru vs nonmem
fh4.fig.Children.Subtitle.String='Excitatory';


if false
    fh4i=bz.inter_wave_pct(wrs_mux_meta,'inhibit',true);
    fh4i.fig.Children.Subtitle.String='Inhibitory';
end
%% FIG 4 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
% [sig,pair]=bz.load_sig_sums_conn_file('pair',true);
fcstats=fc.fc_com_reg_wave.stats(wrs_meta,com_map,reg_com_maps,'odor_only',true);
[barmm,barci,barcnt]=fc.fc_com_reg_wave.sums(fcstats,'odor_only',true);
fh=fc.fc_com_reg_wave.plot(barmm,barci,barcnt,'condense_plot',false,'odor_only',true,'omit_reg_wave',false);


if false
    inter_wave_fh=bz.inter_wave_ext_bars(wrs_mux_meta);  % dur, olf vs Isocortex, Straitum and Midbrain
    %skipped for current manuscript
end
% bz.inter_wave_ext_bars()

% separate wave timing cdf
% wave.mix_single_wave_timing

%>>> jump to TCOM section as needed

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

%% chain
chains_fwd=wave.COM_chain(wrs_meta,com_map,'odor_only',true);
chains_rev=wave.COM_chain(wrs_meta,com_map,'odor_only',true,'reverse',true);

nnz(cellfun(@(x) numel(unique(x)),chains_fwd.cids)>5)
nnz(cellfun(@(x) numel(unique(x)),chains_rev.cids)>5)

wave.chain_stats(chains_fwd,chains_rev,su_meta,wrs_meta);
wave.chain_stats_regs(chains_fwd,su_meta,"len_thresh",6,"odor_only",true)
if false
    out=wave.chain_tag(chains_fwd,'skip_save',true,'len_thresh',6); % per-spk association
    save(fullfile('bzdata','chain_tag_fdr_6.mat'),'out','blame')
else
    load(fullfile('bzdata','chain_tag_fdr_6.mat'),'out')
end
wave.motif_dynamic.single_spike_chains(out)

wave.chains_loops_sc
wave.chain_SC %plot

%% loops
% sums_all
bz.rings.ring_wave_freq(wrs_meta,'denovo',true,'burst',false,'repeats',3); 

load(fullfile('bzdata','sums_ring_stats_all.mat'),'sums_all');
if false
    pstats=bz.rings.rings_time_constant.stats(sums_all,wrs_meta,'load_file',false,'skip_save',true);
else
    load(fullfile('bzdata','loops_stats_fdr_6.mat'),'pstats')
end
wave.motif_dynamic.single_spike_loops(pstats)
bz.rings.loop_occurrence_per_reg_su(sums_all,su_meta,wrs_meta);

if false
    bz.rings.rings_su_wave_tcom_corr(sums_all)
end

%% chained loops

wave.module_motif_asso_composite
wave.chain_loop_stats



%% w/ bursts

wave.chain_sust_tag(chains_fwd,'burstInterval',150)
wave.chain_sust_tag(chains_fwd,'burstInterval',300)
wave.chain_sust_tag(chains_fwd,'burstInterval',600)

% chains, inconsistent (reverse) direction
wave.chain_tag(chains_rev,'rev',true) % per-spk association

rev_out_150=wave.chain_sust_tag(chains_rev,'burstInterval',150,'rev',true);



%% exports

