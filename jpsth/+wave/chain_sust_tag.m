%
% load('chains_mix.mat','chains_uf');
% chains=chains_uf;
% clear chains_uf;

function out=chain_sust_tag(chains,opt)
arguments
    chains
    opt.ccg (1,1) logical = true
end

% for devp

tsidin={[0:150:400,1500:150:1900],[300:150:700,1800:150:2200],[600:150:1000,2100:150:2500]};
% tsidin=evalin('base','in');
out=chain_recur(tsidin,[],[],1);
return

%% DEBUG
% out=chain_alt(chains);
%%
if opt.ccg
    load('sums_conn_10.mat','sums_conn_str');
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
                %                 ts_id=sortrows(ts_id,1);
                tsidin=cell(0);
                for kk=1:numel(cids)
                    tsidin{kk}=ts_id(ts_id(:,3)==kk,1);
                end
                ts=chain_recur(tsidin,[],[],[]);
                if ~isempty(ts)

                    outkey="s"+sessid+"c"+cc;
                    out.("d"+duration).(wid).(outkey).ts=ts;
                    out.("d"+duration).(wid).(outkey).meta={cids,chains.tcoms(cc)};
                    out.("d"+duration).(wid).(outkey).ts_id=ts_id;

                    % TODO: ccg
                    if opt.ccg
                        cursess=chains.sess(cc);
                        sesspath=ephys.sessid2path(cursess);
                        strippath=regexp(sesspath,'(?<=\\).*','match');
                        sesssel=find(contains({sums_conn_str.folder},strippath));
                        ccg=sums_conn_str(sesssel).ccg_sc;
                        ccgid=sums_conn_str(sesssel).sig_con;
                        chainccg=[];
                        for ii=1:numel(cids)-1
                            sigsel=ccgid(:,1)==cids(ii) & ccgid(:,2)==cids(ii+1);
                            if nnz(sigsel)~=1
                                keyboard()
                            end
                            chainccg=[chainccg;ccg(sigsel,:)];
                        end
                        out.("d"+duration).(wid).(outkey).ccgs=chainccg;
                    end
                end
            end
        end
    end
end

for dur=reshape(fieldnames(out),1,[])
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            ccgs=out.(dur{1}).(wv{1}).(lp{1}).ccgs;
            ccg_qual=[];
            for ii=1:size(ccgs,1)
                ccg_qual=[ccg_qual;bz.good_ccg(ccgs(ii,:))];
                %1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising
                %edge 6:falling edge
            end
            out.(dur{1}).(wv{1}).(lp{1}).ccg_qual=ccg_qual;
        end
    end
end

blame=vcs.blame();
save('chain_tag.mat','out','blame')
stats(out);
end

function stats(out)
for dur=reshape(fieldnames(out),1,[])
    perchaindur=struct();
    [perchaindur.size,perchaindur.dur]=deal([]);
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            perchaindur.size=[perchaindur.size;size(out.(dur{1}).(wv{1}).(lp{1}).ts,2)];
            perchaindur.dur=[perchaindur.dur;{diff(out.(dur{1}).(wv{1}).(lp{1}).ts(:,[1,end]),1,2)./30}];
        end
    end
    statss.("d"+dur)=perchaindur;
end

disp([max(cell2mat(statss.dd3.dur)),max(cell2mat(statss.dd6.dur))])

end


function out=chain_recur(in,chainDepth,currIdx,recDepth)
disp("recDepth"+num2str(recDepth))

out=cell(0); %{[depth,id]}
clen=size(in,2);

if isempty(chainDepth)
    chainDepth=1;
    currIdx=1;
end

while currIdx<=numel(in{chainDepth})
    % no branch
    if currIdx<=numel(in{chainDepth})-1 && in{chainDepth}(currIdx+1)-in{chainDepth}(currIdx)<300
        link=[chainDepth,currIdx,in{chainDepth}(currIdx)];
        fprintf('%d%dS>',chainDepth,currIdx)
        rec=chain_recur(in,chainDepth,currIdx+1,recDepth+1);
        fprintf('%d%dS<',chainDepth,currIdx)
%         if chainDepth==1 && currIdx==1
%             keyboard()
%         end
        for ii=1:numel(rec)
            out(end+1)={[link;rec{ii}]};
        end
    else % no branch last ts  %TODO varify time window
        if chainDepth==clen
            out={[chainDepth,currIdx,in{chainDepth}(currIdx)]};
        end
    end
    % branch
    if chainDepth<clen
        nxtd=chainDepth+1;
        nidx=1;
        while nidx<=numel(in{nxtd}) && in{nxtd}(nidx)<in{chainDepth}(currIdx)+24
            nidx=nidx+1;
        end
        if in{nxtd}(nidx)<in{chainDepth}(currIdx)+300 % in window
            link=[chainDepth,currIdx,in{chainDepth}(currIdx)];
            fprintf('%d%dB>',chainDepth,currIdx)
            rec=chain_recur(in,nxtd,nidx,recDepth+1);
            fprintf('%d%dB<',chainDepth,currIdx)
            for ii=1:numel(rec)
                out(end+1)={[link;rec{ii}]};
            end
        end
    end
    if recDepth==1  % TODO: if branch, 13B loops, if not 11S skips
        while currIdx<=numel(in{1})-1 && in{1}(currIdx+1)<in{1}(currIdx)+300
            currIdx=currIdx+1;
        end
        currIdx=currIdx+1;
    else
        break
    end
end
end



