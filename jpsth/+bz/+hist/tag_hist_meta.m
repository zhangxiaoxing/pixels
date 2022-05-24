function [reg,mem_type]=tag_hist_meta(per_type_stats,opt)
arguments
    per_type_stats (1,1) struct
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
persistent meta type criteria
if isempty(meta) || ~strcmp(type,opt.type) || ~strcmp(criteria,opt.criteria)
    meta=ephys.util.load_meta('type',opt.type,'criteria',opt.criteria);
    meta.sessid=int32(cellfun(@(x) ephys.path2sessid(x,'type',opt.type,'criteria',opt.criteria),meta.allpath));
    type=opt.type;
    criteria=opt.criteria;
end
reg_map=containers.Map('KeyType','int32','ValueType','any');
mem_type_map=containers.Map('KeyType','int32','ValueType','any');
for i=1:numel(meta.sessid)
    key=bitshift(meta.sessid(i),16)+int32(meta.allcid(i));
    reg_map(key)=meta.reg_tree(:,i);
    mem_type_map(key)=meta.mem_type(i);
end
reg_c=arrayfun(@(y) arrayfun(@(x) ...
    reg_map(bitshift(int32(per_type_stats.sess(y)),16)+int32(per_type_stats.sess_suids(y,x))),...
    1:2,'UniformOutput',false), (1:numel(per_type_stats.sess))', 'UniformOutput',false);

reg=[cellfun(@(x) x{1}.',reg_c,'UniformOutput',false),cellfun(@(x) x{2}.',reg_c,'UniformOutput',false)];

mem_type_c=arrayfun(@(y) arrayfun(@(x) ...
    mem_type_map(bitshift(int32(per_type_stats.sess(y)),16)+int32(per_type_stats.sess_suids(y,x))),...
    1:2,'UniformOutput',false), (1:numel(per_type_stats.sess))', 'UniformOutput',false);

mem_type=[cellfun(@(x) x{1},mem_type_c),cellfun(@(x) x{2},mem_type_c)];

end

