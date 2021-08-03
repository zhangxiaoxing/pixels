function span=rings_span(opt)
arguments
    opt.ring_size (1,1) double {mustBeMember(opt.ring_size,3:5)}=3
    opt.to_plot (1,1) logical = false
end
persistent meta
if isempty(meta)
    meta=ephys.util.load_meta();
    meta.sess=cellfun(@(x) ephys.path2sessid(x),meta.allpath);
end
load(fullfile('bzdata','rings_bz.mat'),'rings');
rsidx=opt.ring_size-2;
[~,~,ratiomap]=ref.get_pv_sst();
span=[];
for fi=1:size(rings,1)
    if isempty(rings{fi,rsidx})
        continue
    end
    sess=fi;
    sesssel=meta.sess==sess;
    cids=meta.allcid(sesssel);
    regs=meta.reg_tree(5,sesssel);
    
    for ri=1:size(rings{fi,rsidx},1)
        ring_reg=arrayfun(@(x) regs(cids==x),rings{fi,rsidx}(ri,:),'UniformOutput',false);
        if ~all(cellfun(@(x) ratiomap.isKey(x),ring_reg),'all')
            continue
        end
        pv_ratio=cellfun(@(x) ratiomap(char(x)),ring_reg);
        span=[span;sess,rings{fi,rsidx}(ri,:),max(pv_ratio)-min(pv_ratio)];
    end
end
if opt.to_plot
    figure()
    histogram(span(:,end));
end