function [grey_regs,tbl]=getGreyRegs(opt)
arguments
    opt.mincount (1,1) double = 100
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    opt.table (1,1) logical = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end

meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria",opt.criteria,"load_file",false,"skip_stats",true);
BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
CTXsel=strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),'');
switch opt.range
    case 'grey'
        grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));
    case 'CH'
        grey_regs=unique(meta.reg_tree(5,CHsel));
    case 'CTX'
        grey_regs=unique(meta.reg_tree(5,CTXsel));
end
cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
if opt.table
    tbl.regs=grey_regs;
    tbl.cnt=cnt;
end
grey_regs=grey_regs(cnt>opt.mincount);

end