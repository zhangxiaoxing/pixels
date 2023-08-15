
bz.rings.ring_wave_freq(wrs_mux_meta); 
load(fullfile('bzdata','sums_ring_stats_all.mat'),'sums_all');% 1X3
pstats=bz.rings.rings_time_constant.stats(sums_all,wrs_mux_meta,'load_file',false,'skip_save',true);

if false %denovo
    [sschain.out,unfound]=wave.chain_tag.tag(chains_uf,len_thresh,'skip_save',true,'odor_only',true,'extend_trial',false); % per-spk association
else %load
    load(fullfile('binary','chain_tag_all_trl.mat'))
end
sschain.out=out;

disconnected=wave.module_motif_asso_composite(sschain,pstats);
run_length=wave.chain_loop_stats(sschain,pstats,disconnected);

wave.chain_loop_stats;

load(fullfile('binary/','sschain_part.mat'))
sschainp.out=sschain_part;
disconnected=wave.module_motif_asso_composite(sschainp,pstats,trials_dict);
run_length=wave.chain_loop_stats(sschainp,pstats,disconnected);
