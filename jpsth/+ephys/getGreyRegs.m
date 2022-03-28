function grey_regs_=getGreyRegs(opt)
arguments
    opt.mincount (1,1) double = 100
end
persistent grey_regs opt_
if isempty(grey_regs) || ~isequaln(opt_,opt)
    meta=ephys.util.load_meta();
    BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
    CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
    grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));
    cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
    grey_regs=grey_regs(cnt>opt.mincount);
end
grey_regs_=grey_regs;
opt_=opt;
end