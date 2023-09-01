function gen_chain_tag_file(opt)
arguments
    opt.poolsize=2
end
load(fullfile('binary','su_meta.mat'))
fstr=load(fullfile('binary','wrs_mux_meta.mat'));
sel_meta=fstr.wrs_mux_meta;
clear fstr
global_init;
reg_com_maps=wave.get_reg_com_maps(sel_meta);
chains_uf_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'cross_only',false);
wave.chain_tag.tag('poolsize',opt.poolsize,...
    chains_uf_all,3,'skip_save',false,'odor_only',true,...
    'extend_trial',true,'skip_ts_id',true,...
    'filename','chain_tag_all_trl.mat'); % per-spk association
end