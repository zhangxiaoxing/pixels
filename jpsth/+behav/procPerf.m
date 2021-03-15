function out=procPerf(facSeq, opt)
arguments
    facSeq (:,8) double
    opt.mode   (1,:) char {mustBeMember(opt.mode,{'correct','all','error'})} = 'correct'
end
if strcmp(opt.mode, 'error')
    if length(facSeq)>=40
        errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
        out=facSeq(errorsel,:);
    else
        out=[];
    end
else
    if length(facSeq)>=40
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            good=xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0);
            facSeq(i-39:i,10)=good;   % col9->well-train, col10->correct response
            if nnz(good)>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        if strcmp(opt.mode,'correct')
            out=facSeq(all(facSeq(:,9:10),2),:);
        else
            out=facSeq;
        end
    else
        out=[];
    end
end
end