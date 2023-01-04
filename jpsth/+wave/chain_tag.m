%
% WIP
% load('chains_mix.mat','chains_uf');
% chains=chains_uf;
% clear chains_uf;


function out=chain_tag(chains)
arguments
    chains
end

tags=zeros(size(in,1),1);
chain_idx=1;
curr_pre_ptr=1;
curr_chain=[];
starts=[];
ends=[];
spk_cnt=[];
durs=[];
tsize=size(in,1);

while curr_pre_ptr<tsize
    if rem(curr_pre_ptr,50000)==0, fprintf('%06d.\n',curr_pre_ptr);end
    %matching time window, assuming 30kHz
    %assume max 200hz neuron FR, 2spikes × 5 su in 10ms

    if isempty(curr_chain)
        cyc_number_next=rem(in(curr_pre_ptr,2)+1,rsize);
        if cyc_number_next==0, cyc_number_next=rsize;end
        nxtstep=min(curr_pre_ptr+20,tsize);
        syn_win_ubound=find(in((curr_pre_ptr+1):nxtstep,1)>in(curr_pre_ptr,1)+300,1); %first outside window
        if isempty(syn_win_ubound), break;end %TODO use max available instead
        syn_win_lbound=find(in((curr_pre_ptr+1):(syn_win_ubound+curr_pre_ptr-1),1)>in(curr_pre_ptr,1)+24,1); 
    else
        cyc_number_next=rem(in(curr_chain(end),2)+1,rsize);
        if cyc_number_next==0, cyc_number_next=rsize;end
        nxtstep=min(curr_chain(end)+20,tsize);
        syn_win_ubound=find(in((curr_pre_ptr+1):nxtstep,1)>in(curr_chain(end),1)+300,1); %first outside window
        if isempty(syn_win_ubound), break;end %TODO use max available instead
        syn_win_lbound=find(in((curr_pre_ptr+1):(syn_win_ubound+curr_pre_ptr-1),1)>in(curr_chain(end),1)+24,1); 
    end
    if isempty(syn_win_lbound)
        curr_chain=[];
        curr_pre_ptr=curr_pre_ptr+1;
        continue;
    end %matching time window, assuming 30kHz
    diff_post_ptr=find(in((curr_pre_ptr+syn_win_lbound):(curr_pre_ptr+syn_win_ubound-1),2)==cyc_number_next,1); %post unit
    if ~isempty(diff_post_ptr)
        %TODO temp list chain spk
        if isempty(curr_chain), curr_chain=curr_pre_ptr;end
        curr_pre_ptr=curr_pre_ptr+syn_win_lbound-1+diff_post_ptr;
        curr_chain=vertcat(curr_chain,curr_pre_ptr);
    else
        if numel(curr_chain)>rsize
            tags(curr_chain)=chain_idx;
            chain_idx=chain_idx+1;
            starts(end+1)=in(curr_chain(1),2);
            ends(end+1)=in(curr_chain(end),2);
            spk_cnt(end+1)=numel(curr_chain);
            durs(end+1)=in(curr_chain(end),1)-in(curr_chain(1),1);
        end
        curr_chain=[];
        curr_pre_ptr=curr_pre_ptr+1;
%         continue
    end
end

out=struct();
out.tags=sparse(tags);
out.starts=starts;
out.ends=ends;
out.spk_cnt=spk_cnt;
out.durs=durs;
end