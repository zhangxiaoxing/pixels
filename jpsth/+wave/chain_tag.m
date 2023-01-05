%
% WIP
% load('chains_mix.mat','chains_uf');
% chains=chains_uf;
% clear chains_uf;


function out=chain_tag(chains)
arguments
    chains
end
%% build up 
% all_chains=fieldnames(pstats.congru);
waveids=reshape(unique(chains.wave),1,[]);
sesses=reshape(unique(chains.sess),1,[]);

for sessid=sesses
%         sess_chain=all_chains(sess==sessid);
        [spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessid);
    for wid=waveids
        for duration=[3 6]
            if contains(wid,'s1')
                trial_sel=find(trials(:,5)==4 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
            elseif contains(wid,'s2')
                trial_sel=find(trials(:,5)==8 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
            elseif contains(wid,'dur')
                trial_sel=find(trials(:,8)==duration & all(trials(:,9:10)>0,2));
            else
                keyboard();
            end

            for rr=reshape(sess_chain,1,[])
                if ~all(ismember(pstats.congru.(rr{1}).rstats{4},waveid),2) ...
                        || (duration==3 && any(ismember(pstats.congru.(rr{1}).rstats{4},[2,4,8]),"all")) ...
                        || (duration==6 && any(ismember(pstats.congru.(rr{1}).rstats{4},[1,3,7]),"all"))
                    continue
                end
                ts_id=pstats.congru.(rr{1}).ts_id;
                rcids=pstats.congru.(rr{1}).rstats{3};
                for su=rcids
                    sutag="s"+sessid+"w"+wid+"u"+num2str(su);
                    if ~isfield(single_su_multi_chain.("d"+num2str(duration)),sutag)
                        single_su_multi_chain.("d"+num2str(duration)).(sutag)=ts_id(ts_id(:,2)==su & ismember(ts_id(:,5),trial_sel),4:5);
                        ssmr_meta.("d"+num2str(duration)).(sutag)={pstats.congru.(rr{1}).rstats};
                    end
                    single_su_multi_chain.("d"+num2str(duration)).(sutag)(:,end+1)=ts_id(ts_id(:,2)==su & ismember(ts_id(:,5),trial_sel) ,6);
                    ssmr_meta.("d"+num2str(duration)).(sutag)(end+1)={pstats.congru.(rr{1}).rstats};
                end
            end
        end
    end
end














end
function chain_one(in) 
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
