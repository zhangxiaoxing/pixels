%%

global_init;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();

% wrs_mux_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',true,'merge_mux',true,'extend6s',true);

com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false,'odor_only',true);

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
len_thresh=3;
reg_com_maps=cell2struct({tcom3_maps;tcom6_maps},{'tcom3_maps','tcom6_maps'});
chains_uf_all=wave.COM_chain_reg(su_meta,wrs_mux_meta,reg_com_maps);
chains_uf_rev_all=wave.COM_chain_reg(su_meta,wrs_mux_meta,reg_com_maps,'reverse',true);
chains_nm_all=wave.COM_chain_reg(su_meta,wrs_mux_meta,reg_com_maps,'non_mem',true);

fwd_cross=chains_uf_all.cross_reg;
rev_cross=chains_uf_rev_all.cross_reg;
nm_cross=chains_nm_all.cross_reg;
nm_samp=randsample(numel(chains_nm_all.sess),5000);
for fn=reshape(fieldnames(chains_uf_all),1,[])
    chains_uf.(fn{1})=chains_uf_all.(fn{1})(fwd_cross);
    chains_uf_rev.(fn{1})=chains_uf_rev_all.(fn{1})(rev_cross);
    chains_nm.(fn{1})=chains_nm_all.(fn{1})(nm_cross);
    chains_nm_samp.(fn{1})=chains_nm_all.(fn{1})(nm_samp);
end


chains_uf_within=wave.COM_chain(su_meta,wrs_mux_meta,com_map,'odor_only',true);
chains_uf_rev_within=wave.COM_chain(su_meta,wrs_mux_meta,com_map,'reverse',true,'odor_only',true);
% blame=vcs.blame();
% save(fullfile('bzdata','chains_mix.mat'),'chains_uf','chains_uf_rev','blame')

fwd_within=~chains_uf_within.cross_reg;
rev_within=~chains_uf_rev_within.cross_reg;
for fn=reshape(fieldnames(chains_uf_within),1,[])
    chains_uf.(fn{1})=[chains_uf.(fn{1});chains_uf_within.(fn{1})(fwd_within)];
    chains_uf_rev.(fn{1})=[chains_uf_rev.(fn{1});chains_uf_rev_within.(fn{1})(rev_within)];
end

% TODO: remove within-region due to partial-overlap

%%


%%
if false
    load(fullfile('bzdata','sschain_trl_R.mat'))
else
    tic
    [sschain_trl,unfound]=wave.chain_tag.tag(chains_uf,len_thresh,'skip_save',true,'odor_only',true,'extend_trial',true,'skip_ts_id',true); % per-spk association
    toc
    tic
    [sschain_trl_rev,unfound_rev]=wave.chain_tag.tag(chains_uf_rev,len_thresh,'skip_save',true,'odor_only',true,'extend_trial',true,'skip_ts_id',true); % per-spk association
    toc
    tic
    [sschain_trl_nm,unfound_nm]=wave.chain_tag.tag(chains_nm,len_thresh,'skip_save',true,'odor_only',true,'extend_trial',true,'skip_ts_id',true); % per-spk association
    toc
end

%%

[chain_replay,chain_stats,chain_raw]=wave.replay.stats(sschain_trl,'var_len',false);
[chain_replay_rev,chain_stats_rev,chain_raw_incon]=wave.replay.stats(sschain_trl_rev,'var_len',false);
%% vs control
cyy=[chain_stats(1,:),chain_stats_rev(1,:),...
    chain_stats(5,:),chain_stats_rev(5,:),chain_stats(11,:),chain_stats_rev(11,:),...
    chain_stats(12,:),chain_stats_rev(12,:)];
ggn=[size(chain_stats,2),size(chain_stats_rev,2)];

cgg=[ones(ggn(1),1);2*ones(ggn(2),1);...
    3*ones(ggn(1),1);4*ones(ggn(2),1);...
    5*ones(ggn(1),1);6*ones(ggn(2),1);...
    7*ones(ggn(1),1);8*ones(ggn(2),1)];

cmm=arrayfun(@(x) mean(cyy(cgg==x & isfinite(cyy.'))),1:8);
cci=cell2mat(arrayfun(@(x) bootci(100,@(x) mean(x), cyy(cgg==x & isfinite(cyy.'))),1:8,'UniformOutput',false));


figure()
hold on
bar(cmm.','grouped','FaceColor','none','EdgeColor','k')
errorbar(1:numel(cmm),cmm,cci(1,:)-cmm,cci(2,:)-cmm,'k.');
set(gca(),'XTick',1.5:2:10,'XTickLabel',{'Delay','ITI','Before','After'})
title('chains consis-incon')

%% region
wave.replay.region_replay(chain_replay,'reg',"HIP")