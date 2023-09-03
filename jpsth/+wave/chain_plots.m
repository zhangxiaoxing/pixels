%%

%% overall freq
load(fullfile('binary','motif_replay.mat'))
delayfreq=[loops_sums(1,:),chain_sums(1,:)];
prctile(delayfreq,[25,50,75])

[bhat,bci]=binofit(numel(unique([...
cell2mat(arrayfun(@(x) double(chain_replay.session(x)).*100000+double(chain_replay.meta{x,2}),1:size(chain_replay,1),'UniformOutput',false)),...
cell2mat(arrayfun(@(x) double(ring_replay.session(x)).*100000+double(ring_replay.meta{x,2}),1:size(ring_replay,1),'UniformOutput',false))...
])),nnz(wrs_mux_meta.wave_id>0 & wrs_mux_meta.wave_id<7));

figure()
hold on
bar(bhat,'FaceColor','w')
errorbar(1,bhat,bci(1)-bhat,bci(2)-bhat,'k.')
ylim([0,0.2]);
set(gca,'YTick',0:0.1:0.2,'YTickLabel',0:10:20,'XTick',[])
appendfig('tag','motif neuron / wm neuron, chain_plots.m')
%%

numel(chains_uf_all.sess)./sum([numel(chains_uf_all.sess),numel(chains_uf_rev_all.sess)])
numel(chains_nm_all.sess)./sum([numel(chains_nm_all.sess),numel(chains_nm_rev_all.sess)])

wave.chain_stats(chains_uf_all,chains_uf_rev_all,su_meta,'odor_only',true); % static number count
wave.chain_stats(chains_nm_all,chains_nm_rev_all,su_meta,'odor_only',true);


% will not 
if false
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
    toc % ~17 sec

    %% TODO: how to deal with within-region chains?
    tic
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

    % remove within-region due to partial-overlap
    toc % ~3s
end

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
    [sschain_trl_nm_chains_nm_samp,unfound_nm_chains_nm_samp]=wave.chain_tag.tag(chains_nm_samp,len_thresh,'skip_save',true,'odor_only',true,'extend_trial',true,'skip_ts_id',true); % per-spk association

    % [sschain_trl_nm,unfound_nm]=wave.chain_tag.tag(chains_nm,len_thresh,'skip_save',true,'odor_only',true,'extend_trial',true,'skip_ts_id',true); % per-spk association
    toc

tic
[sschain_trl_nm,unfound_nm]=wave.chain_tag.tag(chains_nm,len_thresh,'skip_save',true,'odor_only',true,'extend_trial',true,'skip_ts_id',true); % per-spk association
toc
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


%% vs non-mem motif spike rate
% complete and self-sufficient if data file present


[chain_replay_nm,chain_stats_nm,chain_raw_incon]=wave.replay.stats(sschain_trl_nm_samp,'var_len',false);
cyy=[chain_stats(1,:),chain_stats_nm(1,:),...
    chain_stats(5,:),chain_stats_nm(2,:),chain_stats(11,:),chain_stats_nm(4,:),...
    chain_stats(12,:),chain_stats_nm(5,:)];

ggn=[size(chain_stats,2),size(chain_stats_nm,2)];

cgg=[ones(ggn(1),1);2*ones(ggn(2),1);...
    3*ones(ggn(1),1);4*ones(ggn(2),1);...
    5*ones(ggn(1),1);6*ones(ggn(2),1);...
    7*ones(ggn(1),1);8*ones(ggn(2),1)];

cmm=arrayfun(@(x) median(cyy(cgg==x & isfinite(cyy.'))),1:8);
cci=cell2mat(arrayfun(@(x) bootci(100,@(x) median(x), cyy(cgg==x & isfinite(cyy.'))),1:8,'UniformOutput',false));


figure()
hold on
bar(cmm.','grouped','FaceColor','none','EdgeColor','k')
errorbar(1:numel(cmm),cmm,cci(1,:)-cmm,cci(2,:)-cmm,'k.');
set(gca(),'XTick',1.5:2:10,'XTickLabel',{'Delay','ITI','Before','After'})
title('chains consis vs nonmem 5000 sample')



%% unused control stats 

yy=[chain_stats(1,:),chain_stats_rev(1,:)];
xx=[zeros(size(chain_stats(1,:))),ones(size(chain_stats_rev(1,:)))];
figure()
swarmchart(xx,yy)
set(gca,'YScale','log','YLim',[0.001,100])

mean(chain_stats(1,:))
mean(chain_stats_rev(1,:))

figure
hold on
fwdh=cdfplot(chain_stats(1,:));
revh=cdfplot(chain_stats_rev(1,:));
fwdh.Color='r';
revh.Color='k';
set(gca,'XScale','log')



