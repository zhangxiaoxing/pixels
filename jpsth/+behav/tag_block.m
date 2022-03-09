function dur_resp=tag_block(trials,opt)
arguments
    trials
    opt.wt (1,1) logical = true
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

    if opt.wt
        dur_resp=dur_resp(trials(:,9)~=0,:);
    end
end