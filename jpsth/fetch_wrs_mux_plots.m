% TODO: switch neuron?
% TODO: FC/Chain both-duration congruent 

%% basic stats %%%%%%%%%%%%%%%%%%%%>>
global_init;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();

trials_dict=behav.get_trials_dict('skip_save',true);

% wrs_mux_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',true,'merge_mux',true,'extend6s',true);

% map_cells: mixed_map,olf_map,dur_map
% TODO: cross_thresh hold
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'odor_only',true);
%
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
tcom3_maps=struct();
tcom6_maps=struct();
grey_regs=ephys.getGreyRegs('range','grey');

[fcom3.odor_only.collection,fcom3.odor_only.com_meta]=wave.per_region_COM(...
    com_map,'sel_type','odor_only','com_field','com3');
ureg=intersect(grey_regs,fcom3.odor_only.collection(:,2));
[~,tcidx]=ismember(ureg,fcom3.odor_only.collection(:,2));
tcom3_maps.odor_only=containers.Map(...
    ureg,num2cell(cellfun(@(x) x/4, fcom3.odor_only.collection(tcidx,1))));

[fcom6.odor_only.collection,fcom6.odor_only.com_meta]=wave.per_region_COM(...
    com_map,'sel_type','odor_only','com_field','com6');
ureg=intersect(grey_regs,fcom6.odor_only.collection(:,2));
[~,tcidx]=ismember(ureg,fcom6.odor_only.collection(:,2));
tcom6_maps.odor_only=containers.Map(...
    ureg,num2cell(cellfun(@(x) x/4, fcom6.odor_only.collection(tcidx,1))));

%%


% test for different proportion 
if false
on=nnz(ismember(wrs_mux_meta.wave_id,5:6));
bn=nnz(ismember(wrs_mux_meta.wave_id,1:4));
dn=nnz(ismember(wrs_mux_meta.wave_id,7:8));
alln=numel(wrs_mux_meta.wave_id);

[~,~,p]=crosstab([zeros(on,1);ones(alln-on,1);...
    zeros(bn,1);ones(alln-bn,1);...
    zeros(dn,1);ones(alln-dn,1)],...
    [ones(alln,1);2*ones(alln,1);3*ones(alln,1)])
end


%% wave & stay
wave_n_stay=nnz(ismember(wrs_mux_meta.wave_id,1:6) & wrs_mux_meta.p_olf(:,3)<0.05 & all(wrs_mux_meta.p_olf6<0.05,2));
olf=nnz(ismember(wrs_mux_meta.wave_id,1:6));

%% show case %%%%%%%%%%%%%%%%%%%%%%%%>
% olf %%%%%%%%%%%%%%%%>
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

%<<<<<<<<<<<<<<<<<<<<<<
%%%>dur %%%%%%%%%%%%%%>>
% idx=find(ismember(wrs_mux_meta.wave_id,7:8) & any(wrs_mux_meta.p_dur<1e-4,2));
idx=[2617];
for ii=reshape(idx,1,[])
    scfh=ephys.sens_dur_SC(ii,su_meta,'skip_raster',false,'skip_fill',true);%
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
    scfh=ephys.sens_dur_SC(ii,su_meta,'skip_raster',false,'skip_fill',true);%
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


%% correct error decoding %%%%%%%%%%>>
% svm on neuron firing rates
if false
    odor4odor=pct.pct_decoding_correct_error(wrs_mux_meta,5:6,'lblidx',5,'n_su',50,'rpt',2);% odor
    % dur4odor=pct.pct_decoding_correct_error(wrs_mux_meta,7:8,'lblidx',5,'n_su',50);% odor
    % mux4odor=pct.pct_decoding_correct_error(wrs_mux_meta,1:4,'lblidx',5,'n_su',50);% odor

    % odor4dur=pct.pct_decoding_correct_error(wrs_mux_meta,5:6,'lblidx',8,'n_su',50);% dur
    % dur4dur=pct.pct_decoding_correct_error(wrs_mux_meta,7:8,'lblidx',8,'n_su',50);% dur
    % mux4dur=pct.pct_decoding_correct_error(wrs_mux_meta,1:4,'lblidx',8,'n_su',50);% dur
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
if false
    ephys.plot_switch_cont_decoding(wrs_mux_meta,"type",'olf')
    ephys.plot_switch_cont_decoding(wrs_mux_meta,"type",'dur')
    ephys.plot_switch_cont_decoding(wrs_mux_meta,"type",'mix')
end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%% Sustained vs transient, AUROC
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

while true
    wave_half_half_fh=wave.plot_wave_half_half(wrs_mux_meta,'minr',0.90);
    if ~isempty(wave_half_half_fh)
        break;
    end
end
stats_half_half_fh=wave.COM_half_half_wrs_mux('iti',true);

%% wave 
fh3=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0.1,0.7],'gauss2d',true,'delay',3);
[fh6,imdata]=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0.1,0.7],'gauss2d',true,'delay',6);
if false
fid=fopen(fullfile('binary','upload','F1Q_R_Mean_FR_wave_sort_by_TCOM.json'),'w');
fprintf(fid,jsonencode(imdata));
fclose(fid);
end

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Proportion, TCOM related

[frac_map_cells,pct_bar_fh]=ephys.pct_reg_bars(su_meta,wrs_mux_meta,'xyscale',{'linear','linear'},'only_odor',true); % only need map_cells for tcom-frac corr


%% wave 
fh3=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0,0.7],'gauss2d',true,'delay',3,'xlim',3);
fh6=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0,0.7],'gauss2d',true,'delay',6);

% fraction model
mixed_TCOM_GLM_fh=wave.connectivity_proportion_GLM(frac_map_cells,gather_config.corr_log_log, ...
    'range','grey','data_type','pct-frac','stats_type','percentile');

% TCOM model 3s delay trials only
mixed_TCOM3_GLM_fh=wave.connectivity_proportion_GLM(tcom3_maps,gather_config.corr_ln_log, ...
    'range','grey','data_type','wrs-mux-TCOM3','stats_type','wrs-mux3');

% TCOM model 6s delay trials only
mixed_TCOM6_GLM_fh=wave.connectivity_proportion_GLM(tcom6_maps,gather_config.corr_ln_log, ...
    'range','grey','data_type','wrs-mux-TCOM6','stats_type','wrs-mux6');


%% TCOM and proportion correlation % 4 panel scatters
if false
pct_tcom_fh3=struct();
pct_tcom_fh6=struct();
for typeidx=1:3
    type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
    conn_reg=subsref(["AON","AON","PIR"],struct(type='()',subs={{typeidx}}));
    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom3.(type).collection(:,2));
    ffrac.collection=...
        [num2cell(cellfun(@(x) x(1),frac_map_cells.(type).values(ureg))),...
        ureg,...
        num2cell(ones(numel(ureg),1)*5),...
        num2cell(cellfun(@(x) x(3),frac_map_cells.(type).values(ureg)))];

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
end

%% Functional coupling

%% TODO: COM_CHAIN
% K:\code\jpsth\+wave\COM_chain_SC.m

if false
    inter_wave_fh=bz.inter_wave_ext_bars(wrs_mux_meta);  % dur, olf vs Isocortex, Straitum and Midbrain
    %skipped for current manuscript
end
% bz.inter_wave_ext_bars()

%%% jump to TCOM section as needed
fh4=bz.inter_wave_pct(wrs_mux_meta,'odor_only',true); %congru vs incongru vs nonmem bar lot
fh4.fig.Children.Subtitle.String='Excitatory';
if false
    fh4i=bz.inter_wave_pct(wrs_mux_meta,'inhibit',true);
    fh4i.fig.Children.Subtitle.String='Inhibitory';
end
bz.conn_prob_spatial_dist(sig,pair);
%% FC_Decoding
bz.fccoding.plot_coding(wrs_mux_meta,'dtype','olf')
bz.fccoding.plot_coding(wrs_mux_meta,'dtype','dur')


bz.fc_conn_screen(com_map,pct_meta,'title_suffix','expanded')

%% loops
% sums_all
% tagged data from 'sums_all' -> ring_call.m -> ring_stats.sh -> bz.rings.rings_freq
bz.rings.ring_wave_freq(wrs_mux_meta); 
load(fullfile('bzdata','sums_ring_stats_all.mat'),'sums_all');% 1X3
pstats=bz.rings.rings_time_constant.stats(sums_all,wrs_mux_meta,'load_file',true,'skip_save',true);

% pstats=bz.rings.rings_time_constant.stats(sums_all,wrs_mux_meta,'load_file',false,'skip_save',true);

% bz.rings.rings_reg_pie(sums_all)% 1X3
% bz.rings.rings_freq

bz.rings.loop_occurrence_per_reg_su(sums_all,su_meta,wrs_mux_meta);
if false
    bz.rings.rings_wave_dynamic(sums_all)
    bz.rings.rings_su_wave_tcom_corr(sums_all)
end
%TODO: assembly time constant olf, both, 3s 6s

[~,rings_tag]=bz.rings.rings_time_constant.stats(sums_all,wrs_mux_meta,'load_file',false,'skip_save',true,'compress',true);
load(fullfile('binary','rings_tag.mat'))
[ring_replay,stats_ring,~]=wave.replay.stats(rmfield(rings_tag,"none"),'var_len',true);
ring_replay_tbl=wave.replay.quickconvert(ring_replay);




%% chain join cross-reg, within-reg
% chain_plots.m

% reg_com_maps=cell2struct({tcom3_maps;tcom6_maps},{'tcom3_maps','tcom6_maps'});
% chains_uf_all=wave.COM_chain_reg(su_meta,wrs_mux_meta,reg_com_maps);
% chains_uf_rev_all=wave.COM_chain_reg(su_meta,wrs_mux_meta,reg_com_maps,'reverse',true);
% % TODO: chains nonmem
% % cross
% fwd_cross=chains_uf_all.cross_reg;
% rev_cross=chains_uf_rev_all.cross_reg;
% for fn=reshape(fieldnames(chains_uf_all),1,[])
%     chains_uf.(fn{1})=chains_uf_all.(fn{1})(fwd_cross);
%     chains_uf_rev.(fn{1})=chains_uf_rev_all.(fn{1})(rev_cross);
% end
% 
% chains_uf_within=wave.COM_chain(su_meta,wrs_mux_meta,com_map,'odor_only',true);
% chains_uf_rev_within=wave.COM_chain(su_meta,wrs_mux_meta,com_map,'reverse',true,'odor_only',true);
% 
% 
% % within, using per-su-tcom
% fwd_within=~chains_uf_within.cross_reg;
% rev_within=~chains_uf_rev_within.cross_reg;
% for fn=reshape(fieldnames(chains_uf_within),1,[])
%     chains_uf.(fn{1})=[chains_uf.(fn{1});chains_uf_within.(fn{1})(fwd_within)];
%     chains_uf_rev.(fn{1})=[chains_uf_rev.(fn{1});chains_uf_rev_within.(fn{1})(rev_within)];
% end
% len_thresh=3;
% 
% if false % unhelpful
%     wave.chain_stats(chains_uf,chains_uf_rev,su_meta,'odor_only',true);
% end

% TODO: update shuffle-related parts
% shuf_chains=wave.COM_chain_shuf(wrs_mux_meta,1:100,'odor_only',true)
% % blame=vcs.blame();save('chains_shuf.mat','shuf_chains','blame')
% [fh3,fhf]=wave.chain_stats_regs(chains_uf,su_meta,len_thresh,"odor_only",true);
% set(fh3.Children(1).Children(1),'YLim',[1e-4,2])
% set(fh3.Children(1).Children(2),'YLim',[5e-6,2])
% set(fhf.Children(1).Children(1),'YLim',[1e-3,10])
% set(fhf.Children(1).Children(2),'YLim',[1e-3,10])

% if false %denovo
%     [sschain.out,unfound]=wave.chain_tag.tag(chains_uf,len_thresh,'skip_save',true,'odor_only',true,'extend_trial',false); % per-spk association
% else %load
%     load(fullfile(gather_config.odpath,'Tempdata','chain_tagR.mat'))
% end
% wave.motif_dynamic.single_spike_chains(sschain.out)
% save(fullfile("bzdata","chain_tag_tbl.mat"),"sschain","unfound","blame")
% load(fullfile("bzdata","chain_tag_tbl.mat"),"sschain")



%% replay figure Jun13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% optional import saved data from tempdata folder
% global_init
% load(fullfile(gather_config.odpath,'Tempdata','rings_tag.mat'))

fstr=load(fullfile('binary','chain_tag_all_trl.mat'),'out');
sschain_trl=fstr.out;
[chain_replay,chain_stats,chain_raw]=wave.replay.stats_tbl(sschain_trl,trials_dict,'var_len',false);


if false % skipped due to joint stats
    [chain_replay,chain_stats,chain_raw]=wave.replay.stats_tbl(sschain_trl(cell2mat(sschain_trl.meta(:,3)),:),trials_dict,'var_len',false);
    [cstr,cmat]=wave.replay.stats_replay_sess({chain_raw},'feat_sel',[1 4 5 11 12]);
    fhb=wave.replay.plot_replay_sess_ci(cmat,...
        {'Delay','Test','ITI','Before session','After session',},...
        'title','chains correct trial','ref_line',true,'median_value',true,'stats','median');
    set(gca,'YScale','linear','YLim',[0,6])
end
load(fullfile('binary','rings_tag.mat'))
[ring_replay,ring_stats,ring_raw]=wave.replay.stats(rmfield(rings_tag,'none'),'var_len',true);
rings_tag_nom=rmfield(rmfield(rings_tag,'d6'),'d3');
[ring_nom_replay,ring_nom_stats,ring_nom_raw]=wave.replay.stats(rings_tag_nom,'var_len',true);

if false % skipped due to joint stats
    [rstr,rmat]=wave.replay.stats_replay_sess({ring_raw},'feat_sel',[1 4 5 11 12]);
    fhring=wave.replay.plot_replay_sess_ci(rmat,...
        {'Delay','Test','ITI','Before session','After session',},...
        'title','loops correct trial','ref_line',true,'median_value',true);
    ylim([0.2,10])
end

[jstr,jmat]=wave.replay.stats_replay_sess({chain_raw,ring_raw},'feat_sel',[1 5 11 12]);
fhj=wave.replay.plot_replay_sess_ci(jmat,...
    {'Delay','ITI','Before session','After session',},...
    'title','Motifs correct trial','ref_line',true,'median_value',true);

set(gca(),ylim=[0,3],yscale="linear")
for ii=6:9, fhj.Children.Children(ii).Position(2)=0.25;end


% correct error chain
chain_corr_err=cell2struct({chain_raw.count([1 3 2 5 7 6 8 10 9],:)+eps;...
    chain_raw.time([1 3 2 5 7 6 8 10 9],:);...
    chain_raw.condition; ...
    chain_raw.tag},{'count';'time';'condition';'tag'});
if false % skipped due to joint stats
    [estr,emat]=wave.replay.stats_replay_sess({chain_corr_err});
    fhb=wave.replay.plot_replay_sess_ci(emat,...
        {'Nonpref','Error','Nonpref','Error','Nonpref','Error'},...
        'title','chains delay-prior-after','median_value',false,'ratio_block',3,'ref_p_value',false);
    set(gca,'Ylim',[0,1.75],'YScale','linear')
    srp=[1,signrank(emat(:,1),emat(:,2)),signrank(emat(:,1),emat(:,3)),...
        1,signrank(emat(:,4),emat(:,5)),signrank(emat(:,4),emat(:,6)),...
        1,signrank(emat(:,7),emat(:,8)),signrank(emat(:,7),emat(:,9))];
    bh=findobj(fhb,'-depth',2,'Type','Bar');
    text(bh.XEndPoints,repmat(1.2,1,6), ...
        num2str(srp([2 3 5 6 8 9]).','%.3f'),'HorizontalAlignment','center','VerticalAlignment','top');
end

% correct error loop
ring_corr_err=cell2struct({ring_raw.count([1 3 2 5 7 6 8 10 9],:);...
    ring_raw.time([1 3 2 5 7 6 8 10 9],:);...
    ring_raw.condition; ...
    ring_raw.tag},{'count';'time';'condition';'tag'});
if false % skipped due to joint stats
    [estr,emat]=wave.replay.stats_replay_sess({ring_corr_err});
    fhb=wave.replay.plot_replay_sess_ci(emat,...
        {'Nonpref','Error','Nonpref','Error','Nonpref','Error'},...
        'title','loops delay-prior-after','median_value',false,'ratio_block',3,'ref_p_value',false);
    set(gca,'Ylim',[0,1.2],'YScale','linear')
    srp=[1,signrank(emat(:,1),emat(:,2)),signrank(emat(:,1),emat(:,3)),...
        1,signrank(emat(:,4),emat(:,5)),signrank(emat(:,4),emat(:,6)),...
        1,signrank(emat(:,7),emat(:,8)),signrank(emat(:,7),emat(:,9))];
    bh=findobj(fhb,'-depth',2,'Type','Bar');
    text(bh.XEndPoints,repmat(1.2,1,6), ...
        num2str(srp([2 3 5 6 8 9]).','%.3f'),'HorizontalAlignment','center','VerticalAlignment','top');
    ylabel('Normalized motif spike frequency')
    ylh=yline(1,'r--')
end

[jestr,jemat]=wave.replay.stats_replay_sess({chain_corr_err,ring_corr_err});
fhb=wave.replay.plot_replay_3panel(jemat(:,[1 3 4 6 7 9]),...
    {'Delay','','ITI before','','ITI after',''},...
    'title','loops delay-prior-after','median_value',false);

set(gca,'Ylim',[0.45,1.55],'YScale','linear')
srp=[1,signrank(jemat(:,1),jemat(:,2)),signrank(jemat(:,1),jemat(:,3)),...
    1,signrank(jemat(:,4),jemat(:,5)),signrank(jemat(:,4),jemat(:,6)),...
    1,signrank(jemat(:,7),jemat(:,8)),signrank(jemat(:,7),jemat(:,9))];
bh=findobj(fhb,'-depth',2,'Type','Bar');
text(bh.XEndPoints,repmat(1.2,1,6), ...
    num2str(srp([2 3 5 6 8 9]).','%.3f'),'HorizontalAlignment','center','VerticalAlignment','top');
ylabel('Normalized motif spike frequency')



% chains, vs control
% [chain_replay_anti,chain_stats_anti,chain_raw_anti]=wave.replay.stats(sschain_trl_anti,'var_len',false);
% [cantistr,antimat]=wave.replay.stats_replay_sess({chain_raw_anti});

[chain_replay_rev,chain_stats_rev,chain_raw_incon]=wave.replay.stats(sschain_trl_rev,'var_len',false);
[cinconstr,revmat]=wave.replay.stats_replay_sess({chain_raw_incon});

if false % not suitable due to unmatched network-size
    fhb=wave.replay.plot_replay_cross_sess({cmat(:,1),antimat(:,1),revmat(:,1),...
        cmat(:,3),antimat(:,8),revmat(:,8),...
        cmat(:,4),antimat(:,5),revmat(:,5),...
        cmat(:,5),antimat(:,11),revmat(:,11),...
        cmat(:,6),antimat(:,12),revmat(:,12)},...
        {'','Delay','','','Prior ITI','','','Later ITI','','','Before Session','','','After session',''},...
        'title','chain-consistent-anti-inconsistent','median_value',true);

    [signrank(cmat(:,1),antimat(:,1)), ...
        signrank(cmat(:,3),antimat(:,8)),  ...
        signrank(cmat(:,4),antimat(:,5)),  ...
        signrank(cmat(:,5),antimat(:,11)), ...
        signrank(cmat(:,6),antimat(:,12)); ...

        ranksum(cmat(:,1),revmat(:,1)),  ...
        ranksum(cmat(:,3),revmat(:,8)),  ...
        ranksum(cmat(:,4),revmat(:,5)),  ...
        ranksum(cmat(:,5),revmat(:,11)),  ...
        ranksum(cmat(:,6),revmat(:,12))]


    %  match session
    joinsess=intersect(fieldnames(cinconstr),fieldnames(cstr));
    [csel,~]=ismember(fieldnames(cstr),joinsess);
    [isel,~]=ismember(fieldnames(cinconstr),joinsess);

    [~,cpos]=ismember(joinsess,fieldnames(cstr));
    [~,ipos]=ismember(joinsess,fieldnames(cinconstr));

    fhb=wave.replay.plot_replay_cross_sess({cmat(cpos,1),antimat(cpos,1),revmat(ipos,1),...
        cmat(cpos,3),antimat(cpos,8),revmat(ipos,8),...
        cmat(cpos,4),antimat(cpos,5),revmat(ipos,5),...
        cmat(cpos,5),antimat(cpos,11),revmat(ipos,11),...
        cmat(cpos,6),antimat(cpos,12),revmat(ipos,12)},...
        {'','Delay','','','Prior ITI','','','Later ITI','','','Before Session','','','After session',''},...
        'title','chains-consistent-anti-inconsistent','median_value',true);
end

% chains, vs control v2
cyy=[chain_stats(1,:),chain_stats_rev(1,:),...
    chain_stats(5,:),chain_stats_rev(5,:),chain_stats(11,:),chain_stats_rev(11,:),...
    chain_stats(12,:),chain_stats_rev(12,:)];
ggn=[size(chain_stats,2),size(chain_stats_rev,2)];

cgg=[ones(ggn(1),1);2*ones(ggn(2),1);...
    3*ones(ggn(1),1);4*ones(ggn(2),1);...
    5*ones(ggn(1),1);6*ones(ggn(2),1);...
    7*ones(ggn(1),1);8*ones(ggn(2),1)];

cmm=arrayfun(@(x) median(cyy(cgg==x & isfinite(cyy.'))),1:8);
cci=cell2mat(arrayfun(@(x) bootci(100,@(x) median(x), cyy(cgg==x & isfinite(cyy.'))),1:8,'UniformOutput',false));
if false
figure()
hold on
bar(cmm.','grouped','FaceColor','none','EdgeColor','k')
errorbar(1:numel(cmm),cmm,cci(1,:)-cmm,cci(2,:)-cmm,'k.');
ylim([0,0.22])
set(gca(),'XTick',1.5:2:10,'XTickLabel',{'Delay','ITI','Before','After'})
title('chains consis-incon') % removed anti-direction

% p=kruskalwallis([chain_stats(1,:),chain_stats_anti(1,:),chain_stats_rev(1,:)],[ones(ggn(1),1);2*ones(ggn(2),1);3*ones(ggn(3),1)],'off')
% p=kruskalwallis([chain_stats(8,:),chain_stats_anti(8,:),chain_stats_rev(8,:)],[ones(ggn(1),1);2*ones(ggn(2),1);3*ones(ggn(3),1)],'off')
% p=kruskalwallis([chain_stats(5,:),chain_stats_anti(5,:),chain_stats_rev(5,:)],[ones(ggn(1),1);2*ones(ggn(2),1);3*ones(ggn(3),1)],'off')
% p=kruskalwallis([chain_stats(11,:),chain_stats_anti(11,:),chain_stats_rev(11,:)],[ones(ggn(1),1);2*ones(ggn(2),1);3*ones(ggn(3),1)],'off')
% p=kruskalwallis([chain_stats(12,:),chain_stats_anti(12,:),chain_stats_rev(12,:)],[ones(ggn(1),1);2*ones(ggn(2),1);3*ones(ggn(3),1)],'off')

end


% TODO: loops vs control v2
ryy=[ring_stats(1,:),ring_nom_stats(1,:),...
    ring_stats(5,:),ring_nom_stats(2,:),ring_stats(11,:),ring_nom_stats(4,:),...
    ring_stats(12,:),ring_nom_stats(5,:)];
ggn=[size(ring_stats,2),size(ring_nom_stats,2)];

rgg=[ones(ggn(1),1);2*ones(ggn(2),1);...
    3*ones(ggn(1),1);4*ones(ggn(2),1);...
    5*ones(ggn(1),1);6*ones(ggn(2),1);...
    7*ones(ggn(1),1);8*ones(ggn(2),1)];

rmm=arrayfun(@(x) median(yy(gg==x & isfinite(yy.'))),unique(gg).');
rci=cell2mat(arrayfun(@(x) bootci(100,@(x) median(x), yy(gg==x & isfinite(yy.'))),unique(gg),'UniformOutput',false).');

if false
figure()
% boxplot(yy,gg,'Colors','k','Symbol','c.')
hold on
bar(mm.','grouped','FaceColor','none','EdgeColor','k')
errorbar(1:numel(mm),mm,rci(1,:)-mm.',rci(2,:)-mm.','k.');
ylim([0,1.2])
set(gca(),'XTick',1.5:2:9.5,'XTickLabel',{'Delay','Later','Before','After'})
title('loops congru-nonmem')
ylabel('Motif spike frequencey (Hz)')

end

 

nexttile
hold on
bh=bar([rmm(1:2);rmm(3:4);rmm(5:6);rmm(7:8)],'grouped','EdgeColor','k');
errorbar([bh.XEndPoints],[bh.YEndPoints],...
    rci(1,[1 3 5 7 2 4 6 8])-rmm([1 3 5 7 2 4 6 8]),...
    rci(2,[1 3 5 7 2 4 6 8])-rmm([1 3 5 7 2 4 6 8]),'k.');

set(gca(),'XTick',1:4,'XTickLabel',{'Delay','ITI','Before','After'})
ylabel('Motif spike frequencey (Hz)')
legend(bh,{'Memory','Nonmemory'},'AutoUpdate','off','Location','northoutside','Orientation','horizontal')

rpp=[ranksum(ring_stats(1,:),ring_nom_stats(1,:)),...
    ranksum(ring_stats(5,:),ring_nom_stats(2,:)),...
    ranksum(ring_stats(11,:),ring_nom_stats(4,:)),...
    ranksum(ring_stats(12,:),ring_nom_stats(5,:))]



if false % no working due to unbalanced nunmber and median frequency
    jyy=[cyy,ryy];
    jgg=[cgg;rgg];

    jmm=arrayfun(@(x) median(jyy(jgg==x & isfinite(jyy.'))),unique(jgg).');
    jmean=arrayfun(@(x) mean(jyy(jgg==x & isfinite(jyy.'))),unique(jgg).');
    jci=cell2mat(arrayfun(@(x) bootci(100,@(x) median(x), jyy(jgg==x & isfinite(jyy.'))),unique(jgg),'UniformOutput',false).');

    figure
    hold on
    bh=bar([jmm(1:2);jmm(3:4);jmm(5:6);jmm(7:8)],'grouped','EdgeColor','k')
    errorbar(1:numel(mm),mm,rci(1,:)-mm.',rci(2,:)-mm.','k.');
    ylim([0,1.2])
    set(gca(),'XTick',1:4,'XTickLabel',{'Delay','ITI','Before','After'})
    title('loops congru-nonmem')
    ylabel('Motif spike frequencey (Hz)')
end





%%
if false
wave.chains_time_constant
wave.chains_loops_sc

wave.chain_SC %plot
wave.chain_sust_tag(chains_uf,'burstInterval',150)
wave.chain_sust_tag(chains_uf,'burstInterval',300)
wave.chain_sust_tag(chains_uf,'burstInterval',600)

% chains, inconsistent (reverse) direction
wave.chain_tag(chains_uf_rev,'rev',true) % per-spk association

rev_out_150=wave.chain_sust_tag(chains_uf_rev,'burstInterval',150,'rev',true);
end



%% chained loops
wave.composites


%% exports


pstats=bz.rings.rings_time_constant.stats(sums_all,wrs_mux_meta,'load_file',true,'skip_save',true);

delay_stats=bz.rings.rings_time_constant.stats(sums_all,wrs_mux_meta,'load_file',false,'skip_save',true,'odor_only',true,'delay',true,'iti',false,'remove_non_motif',true);
iti_stats=bz.rings.rings_time_constant.stats(sums_all,wrs_mux_meta,'load_file',false,'skip_save',true,'odor_only',true,'delay',false,'iti',true,'remove_non_motif',true);

dstats=bz.rings.rings_time_constant.plotStats(delay_stats,'odor_only',true,'skip_save',true,'skip_showcase',true);
istats=bz.rings.rings_time_constant.plotStats(iti_stats,'odor_only',true,'skip_save',true,'skip_showcase',true);

blame=vcs.blame();
save(fullfile('binary/','temp0820_delay_iti_loop_stats.mat'),'delay_stats','iti_stats','dstats','istats','blame')

[run_length,covered]=wave.replay.delay_vs_iti(chain_replay,ring_replay_tbl);
fh=wave.replay.plot_delay_vs_iti(run_length);

% for usess=reshape(unique(run_length.delay(:,1)),1,[])
%     c=struct2cell(dstats);
%     cfn=fieldnames(dstats);
%     csel=endsWith(cfn,"s"+usess);
%     sdstats=cell2struct(c(csel),cfn(csel));
% 
%     c=struct2cell(istats);
%     cfn=fieldnames(istats);
%     csel=endsWith(cfn,"s"+usess);
%     sistats=cell2struct(c(csel),cfn(csel));
% 
%     fh=wave.motif_dynamic.composite_loops(sdstats,sistats);
%     title("S"+usess);
% 
%     rr=cell2struct({run_length.delay(run_length.delay(:,1)==usess,:);...
%         run_length.iti(run_length.iti(:,1)==usess,:);...
%         run_length.iti(run_length.delay(:,1)==usess,:)},...
%         {'delay','iti','pre_post'});
% 
%     fh=wave.replay.plot_delay_vs_iti(rr);
%     title("S"+usess);
% 
%     keyboard()
% end
% 


load(fullfile('binary','rings_tag.mat'))
[ring_replay,ring_stats,~]=wave.replay.stats(rmfield(rings_tag,"none"),'var_len',true);
ring_replay_tbl=wave.replay.quickconvert(ring_replay);

fstr=load(fullfile('binary','chain_tag_all_trl.mat'),'out','trials_dict');
[chain_replay,chain_stats,~]=wave.replay.stats_tbl(fstr.out,fstr.trials_dict,'var_len',false);

[run_length,covered]=wave.replay.delay_vs_iti(chain_replay,ring_replay_tbl);

%% before after task duration
onset=[];
offset=[];
for ii=1:116
    t=trials_dict{ii};
    onset=[onset;t(1,1)./30000];
    session_tick=wave.replay.sessid2length(ii);
    offset=[offset;(session_tick-t(end,2))./3000];
end
%%


%%
load(fullfile('binary','su_meta.mat'));
load(fullfile('binary','wrs_mux_meta.mat'));

load(fullfile('binary','sums_ring_stats_all.mat'));
sums=cell(0,4);
for ss=1:3
    for rr=1:size(sums_all{ss},1)
        sess=sums_all{ss}{rr,1};
        cid=sums_all{ss}{rr,3};
        waveid=wrs_mux_meta.wave_id(su_meta.sess==sess & ismember(su_meta.allcid,cid));
        [relat,sel]=bz.rings.ring_wave_type(waveid,"odor_only",true);
        sums=[sums;{ss,rr,relat,sel}];
    end
end

%%
[sig,~]=bz.load_sig_sums_conn_file('pair',false);
sig=bz.join_fc_waveid(sig,wrs_mux_meta.wave_id);
congrucnt=nnz(pct.su_pairs.get_congru(sig.waveid,'odor_only',true));
incongrucnt=nnz(pct.su_pairs.get_incongru(sig.waveid,'odor_only',true));
nonmemcnt=nnz(pct.su_pairs.get_nonmem(sig.waveid));

%%
ring_tag_fstr=load(fullfile('binary','rings_tag_trl.mat'))
ssloop=ring_tag_fstr.ssloop_trl;
load(fullfile('binary','su_meta.mat'));
load(fullfile('binary','wrs_mux_meta.mat'));
greys=ephys.getGreyRegs('range','grey','mincount',0);
regall=cell(0);
for ii=1:size(ssloop,1)
        sess=ssloop.session(ii);
        cid=ssloop.meta{ii,2};
        regs=su_meta.reg_tree(5,su_meta.sess==sess & ismember(su_meta.allcid,cid));
        if ~all(ismember(regs,greys))
            keyboard()
        end
        regall=[regall;{regs}];
end

