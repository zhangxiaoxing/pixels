function out=procPerf(facSeq, opt)
arguments
    facSeq (:,8) double
    opt.mode   (1,:) char {mustBeMember(opt.mode,{'correct','all'})} = 'correct'
end

if length(facSeq)<40 % interrupted sessions
    out=[];
    return
end
if strcmp(opt.mode, 'error')
    errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0); % borrow DNMS rule, as in SBCs
    out=facSeq(errorsel,:);
else
    facSeq(:,9)=0;
    i=40;
    while i<=length(facSeq)
        good=xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0); % borrow DNMS rule, as in SBCs
        facSeq(i-39:i,10)=good;   % col9->well-train, col10->correct response
        if nnz(good)>=30 % 75pct correct rate
            facSeq(i-39:i,9)=1;
        end
        i=i+1;
    end
    facSeq=behav.tag_block(facSeq,'wt',false);
    if strcmp(opt.mode,'correct') % rtn correct trial only
        out=facSeq(all(facSeq(:,9:10),2),:);
    else % 'all'
        out=facSeq;
    end
end
end