function [pvmap,sstmap,ratiomap]=get_pv_sst(opt)
arguments
    opt.raw_ratio (1,1) logical = false
end

persistent pvmap_ sstmap_ ratiomap_
if isempty(pvmap_) || isempty(sstmap_)|| isempty(ratiomap_)
    disp('Rebuild pv-sst ratio map')
    fn='k:\neupix\track_meta\1-s2.0-S0092867417310693-mmc1.xlsx';
    fopts=detectImportOptions(fn);
    fopts.Sheet='count';
    fopts.DataRange=3;
    tbl=readtable(fn,fopts);
    pvmap=containers.Map('KeyType','char','ValueType','double');
    sstmap=containers.Map('KeyType','char','ValueType','double');
    ratiomap=containers.Map('KeyType','char','ValueType','double');

    for i=1:size(tbl,1)
        key=tbl{i,1}{1};
        if ~isempty(key)
            pvmap(key)=tbl{i,4};
            sstmap(key)=tbl{i,6};
            ratiomap(key)=tbl{i,4}./(tbl{i,6}+tbl{i,4});
        end
    end
    pvmap('ACB')=pvmap('ACBcr')+pvmap('ACBsh');
    sstmap('ACB')=sstmap('ACBcr')+sstmap('ACBsh');
    ratiomap('ACB')=pvmap('ACB')./(pvmap('ACB')+sstmap('ACB'));
    
    pvmap_=pvmap;
    sstmap_=sstmap;
    ratiomap_=ratiomap;
else
    disp('Reuse pv-sst ratio map')
    pvmap=pvmap_;
    sstmap=sstmap_;
    ratiomap=ratiomap_;
end