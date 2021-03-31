function out=net_at_depth(sig,pair,depth,opt)
arguments
    sig (1,1) struct
    pair (1,1) struct
    depth (1,1) double {mustBeMember(depth,1:6)}
    opt.homedir (1,:) char = fullfile('K:','code','per_sec');
    opt.min_su (1,1) double = 10
    opt.subsel (1,:) logical = []
    opt.overwrite (1,1) logical = false
    opt.type (1,:) char {mustBeMember(opt.type,{'all','memory','nonmem'})} = 'all'
end
if opt.overwrite || ~isfile(sprintf('net_at_%d_%s.mat',depth,opt.type))
    reg_tree=deblank(h5read(fullfile(opt.homedir,'transient_6.hdf5'),'/reg_tree'));
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
    if isempty(opt.subsel)
        opt.subsel=true(1,size(reg_tree,2));
    end
    
    [regcounts,greg]=groupcounts(reg_tree(depth,opt.subsel)'); % unique region at tree depth
    rsel=cellfun(@(x) ~isempty(x), greg); % skip untagged
    uniqreg=greg(rsel & regcounts>=opt.min_su);
    % regcounts=regcounts(rsel & regcounts>=opt.min_su);
    
    netedge=cellfun(@(x) int32(idmap.reg2ccfid(x)),nchoosek(uniqreg,2));
    netedge=cat(1,netedge,flip(netedge,ndims(netedge)));
    netedge=cat(1,netedge,repmat(unique(netedge(:)),1,2));
    net_count=nan(size(netedge,1),2); % [pair_count, sig_count, ratio]
    if strcmp(opt.type,'memory')
        sigmsel=sig.mem_type>0;
        pairmsel=pair.mem_type>0;
    elseif strcmp(opt.type,'nonmem')
        sigmsel=sig.mem_type==0;
        pairmsel=pair.mem_type==0;
    else
        sigmsel=sig.mem_type>=0;
        pairmsel=pair.mem_type>=0;
    end
    for i=1:size(netedge,1)
        disp(i);
        sigcount=nnz(netedge(i,1)==sig.reg(:,depth,1)...
            & netedge(i,2)==sig.reg(:,depth,2) & sigmsel);
        paircount=nnz(netedge(i,1)==pair.reg(:,depth,1)...
            & netedge(i,2)==pair.reg(:,depth,2) & pairmsel);
        net_count(i,:)=[paircount,sigcount];
    end
    net_count(:,3)=net_count(:,2)./net_count(:,1); % [pair_count, sig_count, ratio]
    netnode=arrayfun(@(x) idmap.ccfid2reg(x),netedge,'UniformOutput',false);
    out=struct;
    out.netnode=netnode;
    out.netedge=netedge;
    out.net_count=net_count;
    save(sprintf('net_at_%d_%s.mat',depth,opt.type),'netedge','net_count','netnode');
else
    out=load(sprintf('net_at_%d_%s.mat',depth,opt.type));
end
end