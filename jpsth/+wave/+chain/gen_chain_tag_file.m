function gen_chain_tag_file(opt)
arguments
    opt.poolsize=2
    opt.nonmem (1,1) logical = false
    opt.shuf (1,1) logical = false
    opt.shufidx = 1
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.shuf_trl (1,1) logical = false
end
switch opt.criteria
    case 'WT'
        load(fullfile('binary','su_meta.mat'),'su_meta')
        fstr=load(fullfile('binary','wrs_mux_meta.mat'));
        sel_meta=fstr.wrs_mux_meta;
        clear fstr
        nm_fn='chain_tag_nonmem_all_trl.mat';
        mem_fn='chain_tag_all_trl.mat';
        shuf_fn='chain_tag_shuf%d.mat';
        shuf_input='bz_ring_shufs.mat';
        shuf_trl_fn='chain_tag_shuftrl%d.mat';        

    case 'Learning'
        su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Learning","load_file",false,"skip_stats",true);
        sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);
        nm_fn='LN_chain_tag_nonmem_all_trl.mat';
        mem_fn='LN_chain_tag_all_trl.mat';
        shuf_fn='LN_chain_tag_shuf%d.mat';
        shuf_input='LN_bz_shufs.mat';
        shuf_trl_fn='LN_chain_tag_shuftrl%d.mat';

    otherwise
        keyboard()
end

global_init;
reg_com_maps=wave.get_reg_com_maps(sel_meta,'criteria',opt.criteria);

% TODO: multi-mode code path

if opt.nonmem
    chains_nm_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'non_mem',true,'cross_only',false,'criteria',opt.criteria);
    %[sschain_trl_nm,unfound_nm]=
    wave.chain_tag.tag(...
        chains_nm_all,3,'skip_save',false,'odor_only',false,...
        'extend_trial',true,'skip_ts_id',true,...
        'filename',nm_fn,'poolsize',opt.poolsize,'criteria',opt.criteria); % per-spk association
else
    if opt.shuf
        load(fullfile('binary',shuf_input),'shufs');
        chains_shuf_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'shuf',true,'shuf_data',shufs{opt.shufidx},'cross_only',false,'non_mem',false,'criteria',opt.criteria);
        wave.chain_tag.tag(...
            chains_shuf_all,3,'skip_save',false,'odor_only',true,...
            'extend_trial',true,'skip_ts_id',true,...
            'filename',fullfile('shufs',sprintf(shuf_fn,opt.shufidx)),'poolsize',opt.poolsize, ...
            'criteria',opt.criteria); % per-spk association
    elseif opt.shuf_trl
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        chains_uf_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'cross_only',false,'criteria',opt.criteria);
        wave.chain_tag.tag(...
            chains_uf_all,3,'skip_save',false,'odor_only',true,...
            'extend_trial',true,'skip_ts_id',true,...
            'shuf_trl',true,'shufidx',opt.shufidx,...
            'filename',fullfile('shufs',sprintf(shuf_trl_fn,opt.shufidx)),'poolsize',opt.poolsize, ...
            'criteria',opt.criteria); % per-spk association


    else
        chains_uf_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'cross_only',false,'criteria',opt.criteria);
        wave.chain_tag.tag(...
            chains_uf_all,3,'skip_save',false,'odor_only',true,...
            'extend_trial',true,'skip_ts_id',true,...
            'filename',mem_fn,'poolsize',opt.poolsize,'criteria',opt.criteria); % per-spk association
    end
end
end