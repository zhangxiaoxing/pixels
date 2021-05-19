function reg=tag_hist_reg(per_type_stats,opt)
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
for i=1:numel(meta.sessid)
    reg_map(bitshift(meta.sessid(i),16)+int32(meta.allcid(i)))=meta.reg_tree(:,i);
end
 reg=arrayfun(@(y) arrayfun(@(x) ...
    reg_map(bitshift(int32(per_type_stats.sess(y)),16)+int32(per_type_stats.sess_suids(y,x))),...
    1:2,'UniformOutput',false), (1:numel(per_type_stats.sess))', 'UniformOutput',false);
end

