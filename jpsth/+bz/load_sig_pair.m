%TODO merge with reg_conn_bz script
function [sig_,pair_]=load_sig_pair(opt)
arguments
    opt.pair (1,1) logical = false
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.prefix (1,:) char = 'BZWT'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.inhibit (1,1) logical = false
end
persistent sig pair type_ criteria_ prefix_
if isempty(sig) || (opt.pair && isempty(pair)) || ~strcmp(opt.type,type_) || ~strcmp(opt.criteria,criteria_) ~strcmp(opt.prefix,prefix_)
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
    if size(fl,1)<9, warning('Files not found');return; end
    sig=struct(); % for significant connect
    sig.suid=cell(0); % cluster id assigned by kilosort, 2nd+ probe prefixed by probe#
    sig.reg=cell(0); % brain region tree
%     sig.wrsp=cell(0); % Wilcoxon rank summation p value
%     sig.selec=cell(0); % selectivity index
    sig.sess=cell(0);
    sig.mem_type=cell(0); % 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched, see epys.get_mem_type


    if opt.pair, pair=sig; end% for all pairs

%     fields={'suid','reg','wrsp','selec','mem_type'};
    fields={'suid','reg','mem_type','per_bin'};
    for fidx=1:size(fl,1)
        disp(fidx);
        fstr=load(fullfile(fl(fidx).folder,fl(fidx).name));
        for fi=fields
            sig.(fi{1}){fidx}=fstr.sig_meta.(fi{1});
            if opt.pair, pair.(fi{1}){fidx}=fstr.pair_meta.(fi{1}); end
        end
        sig.sess{fidx}=repmat(ephys.path2sessid(fstr.pc_stem,'type',opt.type,'criteria',opt.criteria),size(fstr.sig_meta.suid,1),1);
        if opt.pair, pair.sess{fidx}=repmat(ephys.path2sessid(fstr.pc_stem),size(fstr.pair_meta.suid,1),1); end
    end

    for fi=[fields,{'sess'}]
        sig.(fi{1})=cell2mat(sig.(fi{1})');
        if opt.pair, pair.(fi{1})=cell2mat(pair.(fi{1})'); end
    end
end
sig_=sig;
pair_=pair;
type_=opt.type;
prefix_=opt.prefix;
criteria_=opt.criteria;
end