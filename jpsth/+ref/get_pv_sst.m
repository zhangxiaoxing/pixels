function [pvmap,sstmap,ratiomap]=get_pv_sst(opt)
arguments
    opt.raw_ratio (1,1) logical = false
end

persistent pvmap_ sstmap_ ratiomap_
if isempty(pvmap_) || isempty(sstmap_)|| isempty(ratiomap_)
    disp('Rebuild pv-sst ratio map')
    if ispc
        fn='k:\neupix\track_meta\1-s2.0-S0092867417310693-mmc1.xlsx';
    elseif isunix
        fn='~/neupix/track_meta/1-s2.0-S0092867417310693-mmc1.xlsx';
    end
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


%     idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
%     pvsstreg=pvmap.keys();
%     vent_child=cell(0);
%     med_child=cell(0);
%     for reg=reshape(pvsstreg,1,[])
%         if ~idmap.reg2tree.isKey(reg{1})
%             continue
%         end
%         currtree=idmap.reg2tree(reg{1});
%         if find(strcmp(reg{1},currtree))==8
%             if strcmp(currtree{7},'VENT')
%                 vent_child=[vent_child;currtree(8)];
%             elseif strcmp(currtree{7},'MED')
%                 med_child=[med_child;currtree(8)];
%             end
%         end
%     end


%     pvmap('VENT')=sum(cell2mat(pvmap.values(vent_child)));
%     sstmap('VENT')=sum(cell2mat(sstmap.values(vent_child)));
%     ratiomap('VENT')=pvmap('VENT')./(pvmap('VENT')+sstmap('VENT'));
% 
%     pvmap('MED')=sum(cell2mat(pvmap.values(med_child)));
%     sstmap('MED')=sum(cell2mat(sstmap.values(med_child)));
%     ratiomap('MED')=pvmap('MED')./(pvmap('MED')+sstmap('MED'));

    pvmap_=pvmap;
    sstmap_=sstmap;
    ratiomap_=ratiomap;
else
    disp('Reuse pv-sst ratio map')
    pvmap=pvmap_;
    sstmap=sstmap_;
    ratiomap=ratiomap_;
end
