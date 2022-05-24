tbl=readtable(fullfile('K:','neupix','track_meta','structure_tree_safe_2017.csv'));

ccfid2reg=containers.Map('KeyType','uint32','ValueType','any');
ccfid2full=containers.Map('KeyType','uint32','ValueType','any');
reg2ccfid=containers.Map('KeyType','char','ValueType','any');
reg2full=containers.Map('KeyType','char','ValueType','any');
reg2tree=containers.Map('KeyType','char','ValueType','any');
reg2depth=containers.Map('KeyType','char','ValueType','uint32');

for i=1:size(tbl,1)
    ccfid2reg(tbl.id(i))=tbl.acronym(i);
    ccfid2full(tbl.id(i))=tbl.name(i);
    reg2ccfid(tbl.acronym{i})=tbl.id(i);
    reg2full(tbl.acronym{i})=tbl.name(i);
    reg2depth(tbl.acronym{i})=tbl.depth(i);
end

for i=1:size(tbl,1)
    tree=cell(0,0);
    pathstr=strsplit(tbl.structure_id_path{i},'/');
    for j=1:numel(pathstr)
        tree(end+1)=ccfid2reg(str2double(pathstr{j}));
    end
    reg2tree(tbl.acronym{i})=tree;
end


for i=1:size(tbl,1)
    
end


save('reg_ccfid_map.mat','ccfid2full','ccfid2reg','reg2ccfid','reg2full','reg2tree','reg2depth');