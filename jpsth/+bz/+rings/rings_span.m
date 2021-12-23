% TODO return by wave id


function [cross,within]=rings_span(opt)
arguments
    opt.ring_size (1,1) double {mustBeMember(opt.ring_size,3:5)}=3
    opt.to_plot (1,1) logical = false
    opt.memtype (1,:) char {mustBeMember(opt.memtype,{'any','congru','nonmem'})}='any'
    opt.shufid double {mustBeScalarOrEmpty} = []
end

persistent meta rings_shuf
if isempty(meta)
    meta=ephys.util.load_meta();
%     meta.sess=cellfun(@(x) ephys.path2sessid(x),meta.allpath);
end
if ~isempty(opt.shufid) && isempty(rings_shuf)
    load(fullfile('bzdata','rings_bz_shuf.mat'),'rings_shuf');
end

if isempty(opt.shufid)
    load(fullfile('bzdata','rings_bz.mat'),'rings'); % data generated with ring_list_bz.m
else
    rings=rings_shuf{opt.shufid};
end

rsidx=opt.ring_size-2;
[~,~,ratiomap]=ref.get_pv_sst();
if ispc
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
elseif isunix
    idmap=load(fullfile('~','pixels','align','reg_ccfid_map.mat'));
end
[fcom.collection,fcom.com_meta]=wave.per_region_COM('stats_method','mean');
ffrac.collection=ephys.per_region_fraction('memtype','any'); % *100
cross=struct();
cross.meta=[];
cross.reg=cell(0);
within=cross;
cross.pv_ratio=[];
cross.frcom=[];
cross.ffrac=[];
for fi=1:size(rings,1)
    if isempty(rings{fi,rsidx})
        continue
    end
    sess=fi;
    sesssel=meta.sess==sess;
    cids=meta.allcid(sesssel);
    regs=meta.reg_tree(5,sesssel);
    reg_class=meta.reg_tree(2,sesssel);
    memtypes=meta.mem_type(sesssel);
    for ri=1:size(rings{fi,rsidx},1)
        if strcmp(opt.memtype,'nonmem') && ~all(arrayfun(@(x) memtypes(cids==x)==0,rings{fi,rsidx}(ri,:)),'all')
            continue
        elseif strcmp(opt.memtype,'congru') ...
            && ~(all(arrayfun(@(x) ismember(memtypes(cids==x),1:2),rings{fi,rsidx}(ri,:)),'all') ...
            ||all(arrayfun(@(x) ismember(memtypes(cids==x),3:4),rings{fi,rsidx}(ri,:)),'all'))
            continue
        end
        
        ring_class=arrayfun(@(x) reg_class(cids==x),rings{fi,rsidx}(ri,:),'UniformOutput',false);
        if ~all(cellfun(@(x) strcmp(char(x),'CTX'),ring_class),'all')
            continue
        end
        ring_reg=arrayfun(@(x) regs(cids==x),rings{fi,rsidx}(ri,:),'UniformOutput',false);
        if ~all(cellfun(@(x) ratiomap.isKey(x),ring_reg),'all')
            continue
        end
        
        if numel(unique([ring_reg{:}]))==1
            within.meta=[within.meta;sess,rings{fi,rsidx}(ri,:)];
            within.reg=[within.reg;ring_reg{1}];
        else
            pv_ratio=cellfun(@(x) ratiomap(char(x)),ring_reg).*100;
            cross.meta=[cross.meta;sess,rings{fi,rsidx}(ri,:)];
            cross.reg=[cross.reg;ring_reg];
            cross.pv_ratio=[cross.pv_ratio;pv_ratio];
            cross.frcom=[cross.frcom;
                cell2mat(cellfun(@(x) fcom.collection(strcmp(fcom.collection(:,2),x),1),ring_reg)).*0.25];
            cross.ffrac=[cross.ffrac;cell2mat(cellfun(@(x) ffrac.collection(strcmp(ffrac.collection(:,2),x),1),ring_reg)).*100];
        end
    end
end
