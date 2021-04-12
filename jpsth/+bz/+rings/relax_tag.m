function out=relax_tag(in,msize)
% called from \jpsth\+bz\+rings\rings_freq.m

arguments
    in (:,2) double
    msize (1,1) double {mustBeMember(msize,3:5)}
end

tags=zeros(size(in,1),1);
ring_idx=1;
curr_pre_ptr=1;
curr_post_ptr=-1;
curr_ring=[];
starts=[];
ends=[];
spk_cnt=[];
durs=[];
% fprintf('000000');
while curr_pre_ptr<size(in,1)
    if rem(curr_pre_ptr,100)==0, fprintf('%06d.',curr_pre_ptr);end

    cyc_post_pos=rem(in(curr_pre_ptr,2)+1,msize);
    if cyc_post_pos==0, cyc_post_pos=msize;end
    syn_win_ubound=find(in((curr_pre_ptr+1):end,1)>in(curr_pre_ptr,1)+300,1); if isempty(syn_win_ubound), break;end
    syn_win_lbound=find(in((curr_pre_ptr+1):end,1)>in(curr_pre_ptr,1)+15,1); if isempty(syn_win_lbound), break;end %matching time window, assuming 1kHz
    diff_post_ptr=find(in(curr_pre_ptr+1:end,2)==cyc_post_pos,1); %post unit
    if ~isempty(diff_post_ptr) && diff_post_ptr>=syn_win_lbound && diff_post_ptr<syn_win_ubound
        %TODO temp list ring spk
        if isempty(curr_ring), curr_ring=curr_pre_ptr;end
        curr_pre_ptr=curr_pre_ptr+diff_post_ptr;
        curr_ring=vertcat(curr_ring,curr_pre_ptr);
    else
        if numel(curr_ring)>msize
            tags(curr_ring)=ring_idx;
            ring_idx=ring_idx+1;
            starts(end+1)=in(curr_ring(1),2);
            ends(end+1)=in(curr_ring(end),2);
            spk_cnt(end+1)=numel(curr_ring);
            durs(end+1)=in(curr_ring(end),1)-in(curr_ring(1),1);
        end
        curr_ring=[];
        curr_pre_ptr=curr_pre_ptr+1;
%         continue
    end
end

out=struct();
out.tags=tags;
out.starts=starts;
out.ends=ends;
out.spk_cnt=spk_cnt;
out.durs=durs;
end