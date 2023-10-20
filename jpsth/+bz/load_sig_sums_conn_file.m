function [sig_,pair_]=load_sig_sums_conn_file(opt)
arguments
    opt.pair (1,1) logical = false
    opt.fn (1,:) char ...
        {mustBeMember(opt.fn,{'sums_conn_20win.mat','sums_conn_10.mat'})}
    opt.inhibit (1,1) logical = false
    opt.override (1,1) logical = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
global gather_config
if ~isempty(gather_config)
    if gather_config.fc_win==10
        if opt.inhibit
            opt.fn='sums_conn_inhib_10.mat';
        elseif strcmp(opt.criteria,'Learning')
            opt.fn='sums_conn_learning.mat';
        else
            opt.fn='sums_conn_10.mat';
        end
    elseif gather_config.fc_win==20
        if opt.inhibit
            error("Data not ready");
        else
            opt.fn='sums_conn_20win.mat';
        end
    else
        error('Unsupported FC window')
    end
end
assert(isfield(opt,'fn') && ~isempty(opt.fn),'FC window undefined')

persistent sig pair opt_
if isempty(sig) || isempty(pair) || ~isequaln(opt,opt_) || opt.override
    warning(['using FC file ',opt.fn]);
    meta=ephys.util.load_meta('skip_stats',true,'criteria',opt.criteria,'load_file',false,'save_file',false,'adjust_white_matter',true);
    conn_str=load(fullfile("binary",opt.fn));
    flns=fieldnames(conn_str);
    fln=flns{startsWith(flns,'sums_conn')};
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    idmap.reg2ccfid('')=0;
    meta.ccfid=int32(cell2mat(idmap.reg2ccfid.values(meta.reg_tree)));

    sig=struct();
    [sig.suid,sig.reg,sig.sess]=deal([]);
    pair=sig;

    for ii=1:numel(conn_str.(fln))
        t=conn_str.(fln)(ii);
        t.sig_con=reshape(t.sig_con,[],2);
        sess=ephys.path2sessid(t.folder,'criteria',opt.criteria);
        cid=meta.allcid(meta.sess==sess);
        ccfid=meta.ccfid(:,meta.sess==sess);
        sig.suid=[sig.suid;int32(t.sig_con)];
        sig.sess=[sig.sess;repmat(sess,numel(t.sig_con)./2,1)];
        [~,ididx]=ismember(t.sig_con,cid);
        sig.reg=[sig.reg;cat(3,ccfid(:,ididx(:,1)).',ccfid(:,ididx(:,2)).')];

        if opt.pair
            pair_suid=nchoosek(cid,2);
            pair.suid=[pair.suid;int32(pair_suid);int32(pair_suid(:,[2 1]))];
            pair.sess=[pair.sess;repmat(sess,numel(pair_suid),1)];
            [~,ididx]=ismember(pair_suid,cid);
            reg1d=cat(3,ccfid(:,ididx(:,1)).',ccfid(:,ididx(:,2)).');
            pair.reg=[pair.reg;reg1d;reg1d(:,:,[2,1])];
        end
    end
end

sig_=sig;
pair_=pair;
opt_=opt;