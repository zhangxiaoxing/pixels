function [uniqreg,per_reg_stat]=ratio_tree_depth(reg_tree,depth,opt)
arguments
    reg_tree (6,:) cell
    depth (1,1) double {mustBeMember(depth,1:6)}
    opt.subsel (1,:) logical = []
    opt.min_su (1,1) double = 40
end

if isempty(opt.subsel)
    opt.subsel=true(1,size(reg_tree,2));
end

[sust_sel,trans_sel]=basic_stats.get_selective();
[regcounts,greg]=groupcounts(reg_tree(depth,opt.subsel)'); % unique region at tree depth
rsel=cellfun(@(x) ~isempty(x), greg); % skip untagged
uniqreg=greg(rsel & regcounts>=opt.min_su);
regcounts=regcounts(rsel & regcounts>=opt.min_su);
per_reg_stat=nan(numel(uniqreg),3);
for i=1:numel(uniqreg)
    curr_reg=strcmp(reg_tree(depth,:),uniqreg{i}); % region selection
    per_reg_stat(i,:)=[regcounts(i),...
        nnz(sust_sel & curr_reg & opt.subsel),...
        nnz(trans_sel & curr_reg & opt.subsel)];
end
end