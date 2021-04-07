function reg=tag_hist_reg(per_type_stats)
arguments
    per_type_stats (1,1) struct
end
persistent meta
if isempty(meta)
    meta=ephys.util.load_meta();
    meta.sessid=int32(cellfun(@(x) ephys.path2sessid(x),meta.allpath));
end
reg_map=containers.Map('KeyType','int32','ValueType','any');
for i=1:numel(meta.sessid)
    reg_map(bitshift(meta.sessid(i),16)+int32(meta.allcid(i)))=meta.reg_tree(:,i);
end
reg=arrayfun(@(y) arrayfun(@(x) ...
    reg_map(bitshift(int32(per_type_stats.sess(y)),16)+int32(per_type_stats.sess_suids(y,x))),...
    1:2,'UniformOutput',false), (1:numel(per_type_stats.sess))', 'UniformOutput',false);
end
