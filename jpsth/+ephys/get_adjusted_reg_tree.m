function reg_tree=get_adjusted_reg_tree(opt)
arguments
    opt.fn='regs_default.mat'
end
warning(['region tagging file ',opt.fn]);
freg=load(opt.fn,'acronym_list');
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
idmap.reg2tree('')={'','','','','','',''};
tree=idmap.reg2tree.values(freg.acronym_list);
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