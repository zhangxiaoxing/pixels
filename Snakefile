# For quick-start purpose, not fully tested
# Use {repo-root}/jpsth/Snakefile for better compatibility


workdir: "./jpsth"

rule all:
	input:
                "binary/FC_SU_decode_comparison.fig",
                "binary/FC_congru_incon_nonmem.fig",
                "binary/FC_spike_proportion.fig",
                "binary/SC_consistent_inconsistent.fig",
                "binary/SC_rate_vs_distance.fig",
                "binary/behav_previous_trls_glm.fig",
                "binary/chain_freq_w_without_HIP.fig",
                "binary/chain_loop_in_nests.fig",
                "binary/chains_consist_incon.fig",
                "binary/chains_consist_incon_nonmem.fig",
                "binary/chains_occur_per_neuron_per_reg.fig",
                "binary/corr_err_trans_sust_AUC_olf.fig",
                "binary/delay_iti_alt_path_per_spk.fig",
                "binary/delay_iti_motif_spike_per_sec.fig",
                "binary/delay_vs_iti_per_sec.fig",
                "binary/delay_vs_iti_per_sec_w_shuf.fig",
                "binary/fc_coding_congru_incong.fig",
                "binary/frac_neuron_in_motif.fig",
                "binary/hierarchy_time_const_sums.fig",
                "binary/loop_freq_w_without_HIP.fig",
                "binary/loop_occur_vs_shuf.fig",
                "binary/loops_occur_per_neuron_per_reg.fig",
                "binary/motif_freq_chain_mem_vs_nonmem.fig",
                "binary/motif_freq_correct_error.fig",
                "binary/motif_freq_loops_mem_vs_nonmem.fig",
                "binary/motif_freq_prefer_nonpref.fig",
                "binary/motif_spike_total_spike_proportion.fig",
                "binary/nested_loop_time_const_relax_gap.fig",
                "binary/nested_loop_time_constant.fig",
                "binary/nested_motif_replay_after.fig",
                "binary/nested_motif_replay_before.fig",
                "binary/nested_motif_replay_spike_per_sec_w_shuf.fig,",
                "binary/nested_motif_replay_trials.fig",
                "binary/nests_graph_stats.fig",
                "binary/nests_per_region_degree.fig",
                "binary/one_by_one_rmv.fig",
                "binary/sess_motif_spike_per_sec_w_shuf.fig",
                "binary/sschain_time_constant.fig",
                "binary/ssloop_time_constant.fig",
                "binary/su_decoding.fig",
                "binary/su_trans_sustain.fig",
                "binary/chained_loops_extrapolation.fig"

rule trials_dict:
	input:
	output:
		"binary/trials_dict.mat"
	shell:
		'matlab -noFigureWindows -batch "behav.get_trials_dict()"'

rule delay_vs_iti_covered_per_sec:
	input:
		"binary/delay_iti_runlength_covered.mat",
		"binary/trials_dict.mat",
		"binary/delay_iti_runlength_covered_shuf.mat"
	output:
		"binary/delay_vs_iti_per_sec_w_shuf.fig",
		"binary/motif_cover_per_sec.mat"
	shell:
		'matlab -noFigureWindows -batch "wave.replay.delay_vs_iti_per_sec.run_shuf()"'

rule delay_vs_iti_runlength_covered:
	input:
		"binary/motif_replay.mat"
	output:
		"binary/delay_iti_runlength_covered.mat"
	shell:
		'matlab -noFigureWindows -batch "wave.replay.delay_vs_iti()"'

rule delay_vs_iti_runlength_covered_shuf:
	input:
		expand("binary/motif_replay_shuf{sidx}.mat",sidx=list(range(1,101))),
		"binary/trials_dict.mat"
	output:
		"binary/delay_iti_runlength_covered_shuf.mat"
	threads: 52
	shell:
		'matlab -noFigureWindows -batch "wave.replay.delay_vs_iti(shuf=true,poolsize={threads})"'

rule rings_tag_trl:
	input:
		"binary/rings_bz_vs_shuf.mat",
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat"
	output:
		"binary/rings_tag_trl.mat"
	shell:
		'matlab -noFigureWindows -batch "bz.rings.rings_time_constant.stats(load_file=false,skip_save=false,compress=true)"'

rule SC_consistent_inconsistent:
	input:
		"binary/wrs_mux_meta.mat"
	output:
		"binary/SC_consistent_inconsistent.fig"
	shell:
		'matlab -noFigureWindows -batch "fc.fc_com_reg_wave.run_all()"'

rule wrs_mux_meta:
	output:
		"binary/wrs_mux_meta.mat"
	shell:
		'matlab -noFigureWindows -batch "ephys.get_wrs_mux_meta(save_file=true,load_file=false)"'

rule su_meta:
	output:
		"binary/su_meta.mat"
	shell:
		'matlab -noFigureWindows -batch "ephys.util.load_meta(save_file=true,adjust_white_matter=true)"'

rule replay_stats:
	input:
		"binary/trials_dict.mat",
		"binary/chain_tag_all_trl.mat",
		"binary/rings_tag_trl.mat"
	output:
		"binary/motif_replay.mat"
	shell:
		'matlab -noFigureWindows -batch "wave.replay.stats_tbl(skip_save=false)"'

rule gen_all_replay_shuf_stats:
	input:
		expand("binary/motif_replay_shuf{sidx}.mat",sidx=list(range(1,101)))

rule replay_shuf_stats:
	input:
		"binary/trials_dict.mat",
		"binary/chain_tag_shuf{sidx}.mat",
		"binary/ring_tag_shuf{sidx}.mat"
	output:
		"binary/motif_replay_shuf{sidx}.mat"
	shell:
		'matlab -noFigureWindows -batch "wave.replay.stats_tbl(skip_save=false,shuf=true,shufidx={wildcards.sidx})"'

rule delay_iti_motif_per_spike:
	input:
		"binary/motif_replay.mat"
	output:
		"binary/delay_vs_iti_pivot_spike.mat"
	shell:
		'matlab -noFigureWindows -batch "wave.replay.delay_vs_iti_motif_per_spike(skip_save=false)"'

rule delay_iti_motif_spike_plot:
	input:
		"binary/delay_vs_iti_pivot_spike.mat"
	output:
		"binary/delay_iti_alt_path_per_spk.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.replay.plot_pivot()"'

rule motif_nest:
	input:
		"binary/chain_tag_all_trl.mat",
		"binary/rings_tag_trl.mat",
		"binary/trials_dict.mat"
	output:
		"binary/nested_loops_stats.mat",
		"binary/SingleSpikeChainedLoopCid.mat",
		"binary/nests_graph_stats.fig",
		"binary/nests_per_region_degree.fig",
		"binary/chain_loop_in_nests.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.module_motif_asso_composite(skip_save=false)"'

rule nest_loop_runlength:
	input:
		"binary/delay_iti_runlength_covered.mat"
	output:
		"binary/nested_loop_time_constant.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.loop.nested_loop_time_const()"'

rule delay_vs_iti_spike_per_sec:
	input:
		"binary/motif_replay.mat",
		"binary/trials_dict.mat"
	output:
		"binary/delay_iti_motif_spike_per_sec.mat",
		"binary/delay_iti_motif_spike_per_sec.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.replay.delay_vs_iti_motif_spike_per_sec(skip_save=false)"'

rule sschain_time_constant:
	input:
		"binary/motif_replay.mat"
	output:
		"binary/sschain_time_constant.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.chain.sschain_time_const(skip_save=false)"'

rule ssloop_time_constant:
	input:
		"binary/motif_replay.mat"
	output:
		"binary/ssloop_time_constant.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.loop.ssloop_time_const(skip_save=false)"'

rule su_decoding:
	input:
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat"
	output:
		"binary/su_decoding.mat",
		"binary/su_decoding.fig"
	shell:
		'matlab -noFigureWindows -batch "ephys.su_decoding_corr_err(skip_save=false)"'

rule su_trans_sustain_percentage:
	input:
		"binary/wrs_mux_meta.mat"
	output:
		"binary/su_trans_sustain.fig"
	shell:
		'matlab -noFigureWindows -batch "ephys.sust_trans_bar_w_mix(skip_save=false)"'

rule FC_congru_incon_nonmem:
	input:
		"binary/wrs_mux_meta.mat"
	output:
		"binary/FC_congru_incon_nonmem.fig"
	shell:
		'matlab -noFigureWindows -batch "bz.inter_wave_pct(odor_only=true,skip_save=false)"'

rule sust_trans_correct_err_AUC:
	input:
		"binary/wrs_mux_meta.mat"
	output:
		"binary/corr_err_trans_sust_AUC_olf.fig",
		"binary/corr_err_trans_sust_AUC.mat"
	shell:
		'matlab -noFigureWindows -batch "ephys.sust_trans_correct_error(odor_only=true,skip_save=false)"'

rule gen_shuffled_coupling:
	output:
		"binary/bz_ring_shufs.mat"
	threads: 32
	shell:
		'matlab -noFigureWindows -batch "bz.shuffle_conn_bz_alt(poolsize={threads},rpt=100)"'

rule chain_consist_incon:
	input:
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat"
	output:
		"binary/chains_consist_incon.fig",
		"binary/chains_consist_incon_nonmem.fig"
	shell:
		'''
		matlab -noFigureWindows -batch "wave.chain.chains_consist_incon(skip_save=false,shuf=true)"
		matlab -noFigureWindows -batch "wave.chain.chains_consist_incon(skip_save=false,shuf=true,non_mem=true)"
		'''

rule gen_chain_tag_file:
	input:
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat"
	output:
		"binary/chain_tag_all_trl.mat"
	threads: 52
	shell:
		'matlab -noFigureWindows -batch "wave.chain.gen_chain_tag_file(poolsize={threads})"'

rule gen_all_shuf_chain_tag_file:
	input:
		expand("binary/chain_tag_shuf{cidx}.mat",cidx=list(range(1,101)))

rule gen_one_shuf_chain_tag_file:
	input:
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat",
		"binary/bz_ring_shufs.mat"
	output:
		"binary/chain_tag_shuf{cidx}.mat"
	threads: 4
	shell:
		'matlab -noFigureWindows -batch "wave.chain.gen_chain_tag_file(poolsize={threads},shuf=true,shufidx={wildcards.cidx})"'

rule gen_all_shuf_ring_tag_file:
	input:
		expand("binary/ring_tag_shuf{ridx}.mat",ridx=list(range(1,101)))

rule gen_one_shuf_ring_tag_file:
	input:
		"binary/rings_bz_vs_shuf.mat",
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat"
	output:
		"binary/ring_tag_shuf{ridx}.mat"
	shell:
		'matlab -noFigureWindows -batch "bz.rings.rings_time_constant.stats(load_file=false,skip_save=false,compress=true,shuf=true,shufidx={wildcards.ridx})"'

rule motif_freq_w_without_HIP:
	input:
		"binary/motif_replay.mat"
	output:
		"binary/chain_freq_w_without_HIP.fig",
		"binary/loop_freq_w_without_HIP.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.motif_freq_w_without_HIP()"'

rule motif_pref_nonpref_correct_error_freq:
	input:
		"binary/motif_replay.mat"
	output:
		"binary/motif_freq_correct_error.fig",
		"binary/motif_freq_prefer_nonpref.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.motif_prefer_nonprefer()"'
rule gen_one_by_one_rmv_data:
	input:
		"binary/motif_replay.mat"
	output:
		"binary/one_by_one_rmv.mat"
	shell:
		'matlab -noFigureWindows -batch "wave.composite_thin_down.gen_one_by_one_rmv_data()"'

rule plot_one_by_one_rmv_data:
	input:
		"binary/one_by_one_rmv.mat"
	output:
		"binary/one_by_one_rmv.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.composite_thin_down.plot_one_by_one()"'

rule plot_motif_combined_freq_replay:
	input:
		"binary/motif_replay.mat",
		expand("binary/motif_replay_shuf{sidx}.mat",sidx=list(range(1,101))),
		"binary/trials_dict.mat"
	output:
		"binary/nested_motif_replay_spike_per_sec_w_shuf.fig,"
		"binary/sess_motif_spike_per_sec_w_shuf.fig"
	shell:
		'''
		matlab -noFigureWindows -batch "wave.replay.delay_vs_iti_motif_spike_per_sec(skip_save=false,shuf=true,nested=false)"
		matlab -noFigureWindows -batch "wave.replay.delay_vs_iti_motif_spike_per_sec(skip_save=false,shuf=true,nested=true)"
		'''

rule frac_neuron_in_motif:
	input:
		"binary/motif_replay.mat",
		"binary/wrs_mux_meta.mat"
	output:
		"binary/frac_neuron_in_motif.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.frac_neuron_in_motif()"'

rule motif_spike_total_spike_proportion:
	input:
		"binary/delay_vs_iti_pivot_spike.mat"
	output:
		"binary/motif_spike_total_spike_proportion.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.motif_spike_total_spike_proportion()"'

rule gen_shuf_rings:
	input:
		"binary/bz_ring_shufs.mat",
		"binary/sums_conn_10.mat"
	output:
		"binary/rings_bz_vs_shuf.mat"
	threads: 52
	shell:
		'matlab -noFigureWindows -batch "bz.rings.ring_list_bz_alt(poolsize={threads})"'

rule loop_occur_per_su_per_reg:
	input:
		"binary/rings_bz_vs_shuf.mat",
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat",
		"binary/rings_bz_vs_shuf.mat"
	output:
		"binary/loops_occur_per_neuron_per_reg.fig"
	shell:
		'matlab -noFigureWindows -batch "bz.rings.loop_occurrence_per_reg_su()"'

rule chain_occur_per_su_per_reg:
	input:
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat",
		"binary/bz_ring_shufs.mat"
	output:
		"binary/chains_occur_per_neuron_per_reg.fig",
		"binary/chains_occur_per_neuron_per_reg.mat"
	shell:
		'matlab -noFigureWindows -batch "wave.chain.chain_stats_regs()"'

rule gen_fc_coding_data:
	input:
		"binary/sums_conn_10.mat"
	output:
		expand("binary/fccoding/fc_coding_{sess}.mat",sess=list(range(1,117)))
	shell:
		'matlab -noFigureWindows -batch "bz.fccoding.fc_coding_all()"'

rule fc_coding_congru_incon:
	input:
		"binary/wrs_mux_meta.mat",
		expand("binary/fccoding/fc_coding_{sess}.mat",sess=list(range(1,117)))
	output:
		"binary/fc_coding_congru_incong.fig"
	shell:
		'matlab -noFigureWindows -batch "bz.fccoding.plot_coding(odor_only=true,skip_save=false)"'

rule fc_spike_of_all_spike_proportion:
	input:
		"binary/wrs_mux_meta.mat",
		"binary/sums_conn_10.mat"
	output:
		"binary/FC_spike_proportion.mat",
		"binary/FC_spike_proportion.fig"
	shell:
		'matlab -noFigureWindows -batch "bz.FC_spike_proportion()"'

rule fc_su_decode_comparison:
	input:
		"binary/wrs_mux_meta.mat"
	output:
		"binary/FC_SU_decode_comparison.fig"
	shell:
		'matlab -noFigureWindows -batch "bz.fc_su_decode_comparison()"'

rule hierarchy_time_const_sums:
	input:
		"binary/delay_iti_runlength_covered.mat",
		"binary/motif_replay.mat"
	output:
		"binary/hierarchy_time_const_sums.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.rings_wave_dynamic_sums_bar()"'

rule FC_rate_vs_distance:
	input:
		"binary/sums_conn_10.mat"
	output:
		"binary/SC_rate_vs_distance.fig"
	shell:
		'matlab -noFigureWindows -batch "bz.conn_prob_spatial_dist()"'

rule nested_loops_extrapolation:
	input:
		"binary/motif_replay.mat",
		"binary/su_meta.mat"
	output:
		"binary/chainned_loops_thinned.mat",
		"binary/chained_loops_extrapolation.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.chained_loops_extrapol()"'

rule nested_loop_time_const_relax_gap:
	input:
		"binary/delay_iti_runlength_covered.mat"
	output:
		"binary/nested_loop_time_const_relax_gap.fig"
	shell:
		'matlab -noFigureWindows -batch "wave.loop.nested_loop_time_const_relax_gap()"'

rule gen_nonmem_chain_tag_file:
	input:
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat"
	output:
		"binary/chain_tag_nonmem_all_trl.mat"
	threads: 52
	shell:
		'matlab -noFigureWindows -batch "wave.chain.gen_chain_tag_file(poolsize={threads},nonmem=true)"'


rule behav_previous_trls_glm:
	input:
		"binary/trials_dict.mat"
	output:
		"binary/behav_previous_trls_glm.fig"
	shell:
		'matlab -noFigureWindows -batch "behav.previous_trial_glm()"'

rule motif_freq_mem_vs_nonmem:
	input:
		"binary/motif_replay.mat",
		"binary/motif_replay_chain_nonmem.mat"
	output:
		"binary/motif_freq_chain_mem_vs_nonmem.fig",
		"binary/motif_freq_loops_mem_vs_nonmem.fig"
	shell:
		'''
		matlab -noFigureWindows -batch "wave.motif_freq_mem_vs_nonmem(type='chain')"
		matlab -noFigureWindows -batch "wave.motif_freq_mem_vs_nonmem(type='loop')"
		'''


rule nested_motif_showcase:
	input:
		"binary/motif_replay.mat",
		"binary/trials_dict.mat",
		"binary/su_meta.mat",
		"binary/wrs_mux_meta.mat"
	output:
		"binary/nested_motif_replay_trials.fig",
		"binary/nested_motif_replay_before.fig",
		"binary/nested_motif_replay_after.fig"

	shell:
		'matlab -noFigureWindows -batch "wave.replay.motif_seq_activation_replay()"'

rule loop_occurrence_vs_shuf:
	input:
		"binary/wrs_mux_meta.mat"
	output:
		"binary/loop_occur_vs_shuf.fig"
	shell:
		'matlab -noFigureWindows -batch "bz.rings.ring_wave_freq(denovo=true)"'

# /+wave/+loop/motif_seq_activation
# motif_replay_chain_nonmem
