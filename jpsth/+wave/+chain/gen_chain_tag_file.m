function gen_chain_tag_file(opt)
arguments
    opt.poolsize=2
    opt.nonmem (1,1) logical = false
end
load(fullfile('binary','su_meta.mat'))
fstr=load(fullfile('binary','wrs_mux_meta.mat'));
sel_meta=fstr.wrs_mux_meta;
clear fstr
global_init;
reg_com_maps=wave.get_reg_com_maps(sel_meta);
if opt.nonmem
    chains_nm_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'non_mem',true,'cross_only',false);
    %[sschain_trl_nm,unfound_nm]=
    wave.chain_tag.tag(...
        chains_nm_all,3,'skip_save',false,'odor_only',false,...
        'extend_trial',true,'skip_ts_id',true,...
        'filename','chain_tag_nonmem_all_trl.mat','poolsize',opt.poolsize); % per-spk association
else
    chains_uf_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'cross_only',false);
    wave.chain_tag.tag(...
        chains_uf_all,3,'skip_save',false,'odor_only',true,...
        'extend_trial',true,'skip_ts_id',true,...
        'filename','chain_tag_all_trl.mat','poolsize',opt.poolsize); % per-spk association
end
end