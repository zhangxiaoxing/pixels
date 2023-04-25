function [olf_shared,both_shared]=SU_remove_network()
global_init;
load(fullfile('bzdata','chains_mix.mat'),'chains_uf');
% load('chains_shuf.mat','shuf_chains')

chains_uf.uid=arrayfun(@(x) chains_uf.sess(x)*100000+int32(chains_uf.cids{x}), 1:numel(chains_uf.sess),'UniformOutput',false);
len_sel=cellfun(@(x) numel(x),chains_uf.cids)>4;

olf_sel=ismember(chains_uf.wave,["olf_s1","olf_s2"]) & len_sel; % & (chains_uf.reg_sel==curr_tag)
olf_uid=[chains_uf.uid{olf_sel}];
[chain_olf_uid,~]=unique(olf_uid);

both_sel=ismember(chains_uf.wave,["s1d3","s2d3","s1d6","s2d6"]) & len_sel; % & chains_uf.reg_sel==curr_tag
both_uid=[chains_uf.uid{both_sel}];
[chain_both_uid,~]=unique(both_uid);

% dur_sel=ismember(chains_uf.wave,["dur_d3","dur_d6"]) & len_sel; % & chains_uf.reg_sel==curr_tag
% dur_uid=[chains_uf.uid{dur_sel}];
% [chain_dur_uid,~]=unique(dur_uid);

load(fullfile('bzdata','sums_ring_stats_all.mat'),'sums_all');
rstats=bz.rings.rings_reg_pie(sums_all,'plot',false);

lolf_sel=strcmp(rstats(:,9),'olf');
loop_olf_uid=unique([rstats{lolf_sel,10}]);
olf_shared=intersect(chain_olf_uid,loop_olf_uid);

lboth_sel=strcmp(rstats(:,9),'both');
loop_both_uid=unique([rstats{lboth_sel,10}]);
both_shared=intersect(chain_both_uid,loop_both_uid);
end