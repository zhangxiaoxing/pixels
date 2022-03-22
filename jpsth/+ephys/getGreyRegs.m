function grey_regs_=getGreyRegs()
persistent grey_regs
if isempty(grey_regs)
    meta=ephys.util.load_meta();
    BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
    CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
    grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));
    cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
    grey_regs=grey_regs(cnt>100);
end
grey_regs_=grey_regs;
end