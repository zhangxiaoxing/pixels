% TODO: universal wave

global_init;
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false,'odor_only',true);
fh6=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0.1,0.7],'gauss2d',true,'delay',6);