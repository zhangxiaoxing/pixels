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
ch_len=cellfun(@(x) numel(x),chains.cids);

for sessid=sesses
    [spkID,spkTS,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
    for duration=[3 6]
        for wid=waveids
            if (contains(wid,'d3') && duration==6) ...
                || (contains(wid,'d6') && duration ==3)
                continue
            end
            if contains(wid,'s1')
                trial_sel=find(trials(:,5)==4 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
            elseif contains(wid,'s2')
                trial_sel=find(trials(:,5)==8 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
            elseif contains(wid,'dur')
                trial_sel=find(trials(:,8)==duration & all(trials(:,9:10)>0,2));
            else
                keyboard();
            end
            
            sess_indices=reshape(find(chains.sess==sessid & strcmp(chains.wave,wid) & ch_len>4),1,[]);

            for cc=sess_indices
                ts_id=[];
                cids=chains.cids{cc};
%                 disp({sessid,ri});
                for in_chain_pos=1:numel(cids) % TODO, 1:chain_len
                    one_chain_sel=spkID==cids(in_chain_pos);
                    rawts=spkTS(one_chain_sel);

                    ft_sel=strcmp(FT_SPIKE.label,num2str(cids(in_chain_pos)));
                    ft_ts=FT_SPIKE.timestamp{ft_sel};
                    ft_trl_time=FT_SPIKE.time{ft_sel};
                    ft_trl=FT_SPIKE.trial{ft_sel};

                    [~,tspos]=ismember(ft_ts,rawts);
                    ext_time=repmat(-realmax,numel(rawts),1);
                    ext_time(tspos)=ft_trl_time;

                    ext_trl=repmat(-realmax,numel(rawts),1);
                    ext_trl(tspos)=ft_trl;

                    ts_id=cat(1,ts_id,[rawts,... % 1
                        repmat(cids(in_chain_pos),numel(rawts),1),...  % 2
                        ones(numel(rawts),1)*in_chain_pos,...  % 3
                        ext_time,...  % 4
                        ext_trl]); % 5
                end
                % optional remove non-wave spikes
                ts_id=ts_id(ts_id(:,4)>=1 & ts_id(:,4)<(duration+1) & ismember(ts_id(:,5),trial_sel),:);
                %
                ts_id=sortrows(ts_id,1);
                ch_tags=chain_one(ts_id(:,[1 3]),numel(chains.cids{cc}));
                ts_id=[ts_id,ch_tags.tags]; % join TS, chain tag % 6

                % TODO:output
            end
        end
    end
end


end



function out=chain_one(in,chain_len,opt)
arguments
    in
    chain_len
    opt.max_win (1,1) double = 10
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
prevptr=0;
% TODO: trial necessary?
while curr_pre_ptr<tsize
    prevptr=curr_pre_ptr;
    if rem(curr_pre_ptr,5000)==0, fprintf('%06d.\n',curr_pre_ptr);end
    %matching time window, assuming 30kHz
    %assume max 200hz neuron FR, 2spikes × 5 su in 10ms
    if isempty(curr_chain)
        nxtstep=min(curr_pre_ptr+2*chain_len,tsize); % 200hz->2spk/su/10ms
        nxtptr=curr_pre_ptr+find(in(curr_pre_ptr:nxtstep,2)==1,1)-1;
        if isempty(nxtptr)
            curr_pre_ptr=curr_pre_ptr+1;
            continue
        else
            curr_pre_ptr=nxtptr;
        end
        cyc_number_next=2;
        
        syn_win_ubound=find(in((curr_pre_ptr+1):nxtstep,1)>in(curr_pre_ptr,1)+300,1); %first outside window
        if isempty(syn_win_ubound)
            syn_win_ubound=find(in((curr_pre_ptr+1):end,1)>in(curr_pre_ptr,1)+300,1); %first outside window
            if isempty(syn_win_ubound)
                break;
            end
        end %TODO use max available instead
        syn_win_lbound=find(in((curr_pre_ptr+1):(syn_win_ubound+curr_pre_ptr-1),1)>in(curr_pre_ptr,1)+24,1);
    else
        cyc_number_next=in(curr_chain(end),2)+1;
        if cyc_number_next>chain_len
            syn_win_lbound=[];
        else
            nxtstep=min(curr_pre_ptr+2*chain_len,tsize);
            syn_win_ubound=find(in((curr_pre_ptr+1):nxtstep,1)>in(curr_chain(end),1)+300,1); %first outside window
            if isempty(syn_win_ubound)
                syn_win_ubound=find(in((curr_pre_ptr+1):end,1)>in(curr_pre_ptr,1)+300,1); %first outside window
                if isempty(syn_win_ubound)
                    break;
                end
            end %TODO use max available instead
            syn_win_lbound=find(in((curr_pre_ptr+1):(syn_win_ubound+curr_pre_ptr-1),1)>in(curr_chain(end),1)+24,1);
        end
    end
    if isempty(syn_win_lbound)
        curr_chain=[];
        curr_pre_ptr=curr_pre_ptr+1;
        continue;
    end %matching time window, assuming 30kHz
    diff_post_ptr=find(in((curr_pre_ptr+syn_win_lbound):(curr_pre_ptr+syn_win_ubound-1),2)==cyc_number_next,1,'last'); %post unit
    if ~isempty(diff_post_ptr)
        %TODO temp list chain spk
        if isempty(curr_chain), curr_chain=curr_pre_ptr;end
        curr_pre_ptr=curr_pre_ptr+syn_win_lbound-1+diff_post_ptr;
        curr_chain=vertcat(curr_chain,curr_pre_ptr);
    else
        if numel(curr_chain)==chain_len
            tags(curr_chain)=chain_idx;
            chain_idx=chain_idx+1;
        end
        curr_chain=[];
        curr_pre_ptr=curr_pre_ptr+1;
    end
end

out=struct();
out.tags=sparse(tags);
end
