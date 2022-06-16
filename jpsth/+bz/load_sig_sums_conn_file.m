function [sig_,pair_]=load_sig_sums_conn_file(opt)
arguments
    opt.pair (1,1) logical = false
    opt.fn (1,:) char ...
        {mustBeMember(opt.fn,{'sums_conn_20win.mat','sums_conn_10.mat'})}
end
global gather_config
if ~isempty(gather_config)
    if gather_config.fc_win==10
        opt.fn='sums_conn_10.mat';
    elseif gather_config.fc_win==20
        opt.fn='sums_conn_20win.mat';
    else
        error('Unsupported FC window')
    end
end
assert(isfield(opt,'fn') && ~isempty(opt.fn),'FC window undefined')

persistent sig pair opt_
if isempty(sig) || isempty(pair) || ~isequaln(opt,opt_)
    warning(['using FC file ',opt.fn]);
    meta=ephys.util.load_meta();
    conn_str=load(opt.fn);
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    idmap.reg2ccfid('')=0;
    meta.ccfid=int32(cell2mat(idmap.reg2ccfid.values(meta.reg_tree)));

    sig=struct();
    [sig.suid,sig.reg,sig.sess]=deal([]);
    pair=sig;
    for ii=1:numel(conn_str.sums_conn_str)
        t=conn_str.sums_conn_str(ii);
        sess=ephys.path2sessid(t.folder);
        cid=meta.allcid(meta.sess==sess);
        ccfid=meta.ccfid(:,meta.sess==sess);
        sig.suid=[sig.suid;int32(t.sig_con)];
        sig.sess=[sig.sess;repmat(sess,size(t.sig_con,1),1)];
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