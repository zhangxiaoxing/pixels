function reg_tree=get_adjusted_reg_tree(opt)
arguments
    opt.adjust_white_matter (1,1) logical
end

fstr=load(fullfile('binary','regs_default_and_adjust_220608.mat'));

if opt.adjust_white_matter
    acronym_list=fstr.a(:,2);
    warning('Adjusted white matter regions');
else
    acronym_list=fstr.a(:,1);
    warning('Skipped white matter regions');
end
%>>>>>>>>>>>>flatten nested cells>>>>>>>>>>
cellsel=cellfun(@(x) isa(x,'cell'),acronym_list);
acronym_list(cellsel)=cellfun(@(x) x{1},acronym_list(cellsel),'UniformOutput',false);
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
idmap.reg2tree('')={'','','','','','',''};
tree=idmap.reg2tree.values(acronym_list);
reg_tree=cell(6,numel(tree));
for ii=1:numel(tree)
    if numel(tree{ii})>=8
        reg_tree(:,ii)=tree{ii}(3:8);
    else
        reg_tree(:,ii)=repmat({''},1,6);
        reg_tree(1:numel(tree{ii})-2,ii)=tree{ii}(3:end);
    end
end
end