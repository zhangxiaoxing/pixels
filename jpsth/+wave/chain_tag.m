%
% WIP
% load('chains_mix.mat','chains_uf');
% chains=chains_uf;
% clear chains_uf;

function out=chain_tag(chains)
arguments
    chains
end
%% DEBUG
% out=chain_alt(chains);
%%
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
%                 ts_id=sortrows(ts_id,1);
%                 keyboard();
                ts=chain_alt(ts_id(:,[1 3]));
                if ~isempty(ts)
                    outkey="s"+sessid+"c"+cc;
                    out.("d"+duration).(wid).(outkey).ts=ts;
                    out.("d"+duration).(wid).(outkey).meta={chains.cids(cc),chains.tcoms(cc)};
                end
            end
        end
    end
end
stats(out);

end

function stats(out)
for dur=reshape(fieldnames(out),1,[])
    perd=[];
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            perd=[perd;size(out.(dur{1}).(wv{1}).(lp{1}).ts,2),mean(diff(out.(dur{1}).(wv{1}).(lp{1}).ts(:,[1,end]),1,2))./30];
        end
    end
    statss.("d"+dur)=perd;
end
end


function out=chain_alt(in)
out=[];
clen=max(in(:,2));
persu=struct.empty(clen,0);
for ii=1:clen
    persu(ii).tsid=in(in(:,2)==ii,:);
    persu(ii).idx=1;
end

curd=1;
while true
    nxtd=curd+1;
    if persu(nxtd).idx>size(persu(nxtd).tsid,1) || persu(curd).idx>size(persu(curd).tsid,1)
        break
    end
    if persu(nxtd).tsid(persu(nxtd).idx,1)<persu(curd).tsid(persu(curd).idx,1)+24 % ahead of window
        persu(nxtd).idx=persu(nxtd).idx+1;
    elseif persu(nxtd).tsid(persu(nxtd).idx,1)>persu(curd).tsid(persu(curd).idx,1)+300 % beyond window
        persu(curd).idx=persu(curd).idx+1;
        if curd~=1
            curd=curd-1;
        end
    else % in window
        if nxtd==clen
            out=[out;arrayfun(@(x) persu(x).tsid(persu(x).idx,1),1:clen)];
            persu(nxtd).idx=persu(nxtd).idx+1;
        else
            curd=nxtd;
        end
    end
end
end



