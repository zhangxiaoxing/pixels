%TODO merge with reg_conn_bz script
function [sig_,pair_]=load_sig_pair(opt)
arguments
    opt.pair (1,1) logical = false
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.prefix (1,:) char = 'BZWT'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.inhibit (1,1) logical = false
    opt.wf_criteria (1,1) logical = false
    opt.CTX (1,1) logical = true
end
persistent sig pair opt_
if isempty(sig) ...
        || (opt.pair && isempty(pair))...
        || ~isequaln(opt,opt_)
        
    if strcmp(opt.criteria,'Learning') && strcmp(opt.type,'MY')
        fl=dir(fullfile('mydata',sprintf('%s_conn_w_reg_my_learning_*.mat',opt.prefix)));
    elseif strcmp(opt.type,'MY')
        fl=dir(fullfile('mydata',sprintf('%s_conn_w_reg_*.mat',opt.prefix)));
    else
        fl=dir(fullfile('bzdata',sprintf('%s_conn_w_reg_*.mat',opt.prefix)));
        inhibitsel=contains({fl.name},'inhibitory');
        if opt.inhibit
            fl=fl(inhibitsel);
        else
            fl=fl(~inhibitsel);
        end
    end
    
    fprintf('Total sessions %d\n',size(fl,1));
%     if size(fl,1)<9, warning('Files not found');return; end
    sig=struct(); % for significant connect
    sig.suid=cell(0); % cluster id assigned by kilosort, 2nd+ probe prefixed by probe#
    sig.reg=cell(0); % brain region tree
    %     sig.wrsp=cell(0); % Wilcoxon rank summation p value
    %     sig.selec=cell(0); % selectivity index
    sig.sess=cell(0);
    sig.mem_type=cell(0); % 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched, see epys.get_mem_type
    sig.wf_good=cell(0);
    
    if opt.pair, pair=sig; end% for all pairs
    
    %     fields={'suid','reg','wrsp','selec','mem_type'};
    fields={'suid','reg','mem_type','per_bin','wf_good'};
    for fidx=1:size(fl,1)
        disp(fidx);
        fstr=load(fullfile(fl(fidx).folder,fl(fidx).name));
        if opt.wf_criteria
            for fi=fields
                sig.(fi{1}){fidx}=fstr.sig_meta.(fi{1})(all(fstr.sig_meta.wf_good,2),:,:);
                if opt.pair, pair.(fi{1}){fidx}=fstr.pair_meta.(fi{1})(all(fstr.sig_meta.wf_good,2),:,:); end
            end
            sig.sess{fidx}=repmat(ephys.path2sessid(fstr.pc_stem,'type',opt.type,'criteria',opt.criteria),nnz(all(fstr.sig_meta.wf_good,2)),1);
            if opt.pair, pair.sess{fidx}=repmat(ephys.path2sessid(fstr.pc_stem),nnz(all(fstr.sig_meta.wf_good,2)),1); end
            
        else
            for fi=fields
                sig.(fi{1}){fidx}=fstr.sig_meta.(fi{1});
                if opt.pair, pair.(fi{1}){fidx}=fstr.pair_meta.(fi{1}); end
            end
            sig.sess{fidx}=repmat(ephys.path2sessid(fstr.pc_stem,'type',opt.type,'criteria',opt.criteria),size(fstr.sig_meta.suid,1),1);
            if opt.pair, pair.sess{fidx}=repmat(ephys.path2sessid(fstr.pc_stem),size(fstr.pair_meta.suid,1),1); end
        end
    end
    
    for fi=[fields,{'sess'}]
        sig.(fi{1})=cell2mat(sig.(fi{1})');
        if opt.pair, pair.(fi{1})=cell2mat(pair.(fi{1})'); end
    end
    
end
if opt.CTX
    %sig
    ctxsigsel=all(sig.reg(:,2,:)==688,3);
    sig_=struct();
    for f=reshape(fieldnames(sig),1,[])
        fs=f{1};
        if ismatrix(sig.(fs))
            sig_.(fs)=sig.(fs)(ctxsigsel,:);
        else
            sig_.(fs)=sig.(fs)(ctxsigsel,:,:);
        end
    end
    %pair
    if opt.pair
        ctxpairsel=all(pair.reg(:,2,:)==688,3);
        pair_=struct();
        for f=reshape(fieldnames(sig),1,[])
            fs=f{1};
            if ismatrix(sig.(fs))
                pair_.(fs)=pair.(fs)(ctxpairsel,:);
            else
                pair_.(fs)=pair.(fs)(ctxpairsel,:,:);
            end
        end
    else
        pair_=[];    
    end
else
    sig_=sig;
    pair_=pair;
end
opt_=opt;
end