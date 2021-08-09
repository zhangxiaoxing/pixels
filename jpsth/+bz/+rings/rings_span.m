function [cross,within]=rings_span(opt)
arguments
    opt.ring_size (1,1) double {mustBeMember(opt.ring_size,3:5)}=3
    opt.to_plot (1,1) logical = false
    opt.memtype (1,:) char {mustBeMember(opt.memtype,{'any','congru','nonmem'})}='any'
end
persistent meta
if isempty(meta)
    meta=ephys.util.load_meta();
    meta.sess=cellfun(@(x) ephys.path2sessid(x),meta.allpath);
end
load(fullfile('bzdata','rings_bz.mat'),'rings');
rsidx=opt.ring_size-2;
[~,~,ratiomap]=ref.get_pv_sst();
cross=struct();
cross.meta=[];
cross.reg=cell(0);
within=cross;
cross.pv_ratio=[];

for fi=1:size(rings,1)
    if isempty(rings{fi,rsidx})
        continue
    end
    sess=fi;
    sesssel=meta.sess==sess;
    cids=meta.allcid(sesssel);
    regs=meta.reg_tree(5,sesssel);
    reg_class=meta.reg_tree(1,sesssel);
    memtypes=meta.mem_type(sesssel);
    %TODO within region
    for ri=1:size(rings{fi,rsidx},1)
        if strcmp(opt.memtype,'nonmem') && ~all(arrayfun(@(x) memtypes(cids==x)==0,rings{fi,rsidx}(ri,:)),'all')
            continue
        elseif strcmp(opt.memtype,'congru') ...
            && ~(all(arrayfun(@(x) ismember(memtypes(cids==x),1:2),rings{fi,rsidx}(ri,:)),'all') ...
            ||all(arrayfun(@(x) ismember(memtypes(cids==x),3:4),rings{fi,rsidx}(ri,:)),'all'))
            continue
        end
        
        ring_class=arrayfun(@(x) reg_class(cids==x),rings{fi,rsidx}(ri,:),'UniformOutput',false);
        if ~all(cellfun(@(x) strcmp(char(x),'CH'),ring_class),'all')
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
            pv_ratio=cellfun(@(x) ratiomap(char(x)),ring_reg);
            cross.meta=[cross.meta;sess,rings{fi,rsidx}(ri,:)];
            cross.reg=[cross.reg;ring_reg];
            cross.pv_ratio=[cross.pv_ratio;pv_ratio];
        end
    end
end
