function out=procPerf(facSeq, opt)
arguments
    facSeq double
    opt.mode   (1,:) char {mustBeMember(opt.mode,{'correct','all'})} = 'all'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
% assert(~(strcmp(opt.criteria,'WT') && ~size(facSeq,2)==8),"Dimension mismatch");

if length(facSeq)<40 % interrupted sessions
    out=[];
    return
end

facSeq(:,9)=0;
i=40;
while i<=length(facSeq)
    switch opt.criteria
        case 'WT'
            good=xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0); % borrow DNMS rule, as in SBCs
            facSeq(i-39:i,10)=good;   % col9->well-train, col10->correct response
            if nnz(good)>=30 % 75pct correct rate
                facSeq(i-39:i,9)=1;
            end
        case 'Learning'
            good=facSeq(i-39:i,5)~=facSeq(i-39:i,6) & facSeq(i-39:i,7)>0; % borrow DNMS rule, as in SBCs
            if nnz(good)>=15 % 75pct engage
                facSeq(i-39:i,9)=1;
            end
        otherwise
            warning("Unfinished code path")
            keyboard()
    end
    i=i+1;
end
%     block_tag=behav.tag_block(facSeq,'wt',false);
if strcmp(opt.mode,'correct') % rtn correct trial only
    out=facSeq(all(facSeq(:,9:10),2),:);
else % 'all'
    out=facSeq;
end

end