function dur_resp=tag_block(trials,opt)
arguments
    trials
    opt.wt (1,1) logical = true
    opt.block_seq (1,1) logical = false
end

    dur_resp=trials(:,[8,10,8]);
    blk1end=find(diff(dur_resp(:,1))~=0,1);
    revtag=4;
    for i=blk1end:-1:1
        dur_resp(i,3)=revtag;
        revtag=revtag-1;
    end
    
    for i=(blk1end+1):size(dur_resp,1)
        if dur_resp(i,1)~=dur_resp(i-1,1)
            fwd_tag=1;
        else
            fwd_tag=fwd_tag+1;
        end
        dur_resp(i,3)=fwd_tag;
    end

    if opt.block_seq
        bidx=1;
        dur_resp(:,4)=0;
        for jj=1:size(dur_resp,1)
            dur_resp(jj,4)=bidx;
            if dur_resp(jj,3)==4 && dur_resp(jj,1)==6
                bidx=bidx+1;
            end
        end
    end

    if opt.wt
        dur_resp=dur_resp(trials(:,9)~=0,:);
    end
end