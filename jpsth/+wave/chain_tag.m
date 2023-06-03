%
% load(fullfile('bzdata','chains_mix.mat'),'chains_uf');
% chains=chains_uf;
% clear chains_uf;

function [out,notfound]=chain_tag(chains,opt)
arguments
    chains
    opt.ccg (1,1) logical = true
    opt.rev (1,1) logical = false
    opt.shuf_trl (1,1) logical = false
    opt.shuf_idx (1,1) double = 0
    opt.per_reg_wave (1,1) logical = true
    opt.len_thresh (1,1) double = 5
    opt.skip_save (1,1) logical = false
    opt.odor_only (1,1) logical = false
    opt.extend_trial (1,1) logical = false
    opt.skip_ts_id (1,1) logical = false
    opt.DEBUG (1,1) logical = false
    opt.anti_dir (1,1) logical = false
end

if opt.shuf_trl  
    assert(opt.shuf_idx>0,"index is necessary for shuffled data")
end

ch_len=cellfun(@(x) numel(x),chains.cids);
waveids=reshape(unique(chains.wave(ch_len>opt.len_thresh)),1,[]);
sesses=reshape(unique(chains.sess(ch_len>opt.len_thresh)),1,[]);

if opt.ccg
    load('sums_conn_10.mat','sums_conn_str');
end


notfound=cell(0);
%% build chains
% all_chains=fieldnames(pstats.congru);

for sessid=sesses
    if opt.DEBUG && sessid>33
        break
    end
    [~,~,trials,~,~,~]=ephys.getSPKID_TS(sessid,'keep_trial',false);
    for duration=[3 6]
        for wid=waveids
            if (contains(wid,'d3') && duration==6) ...
                || (contains(wid,'d6') && duration ==3)
                continue
            end
            if contains(wid,'s1')
                trial_sel=find(trials(:,5)==4 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
                outid="olf_s1";
            elseif contains(wid,'s2')
                trial_sel=find(trials(:,5)==8 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
                outid="olf_s2";
            elseif contains(wid,'dur')
                if opt.odor_only
                    continue
                else
                    trial_sel=find(trials(:,8)==duration & all(trials(:,9:10)>0,2));
                end
            end
            
            if ~opt.odor_only
                outid=wid;
            end

            sess_indices=reshape(find(chains.sess==sessid & strcmp(chains.wave,wid) & ch_len>=opt.len_thresh),1,[]);
            if isempty(sess_indices), continue;end
            [spkID,spkTS,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
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
                if ~opt.extend_trial
                    ts_id=ts_id(ts_id(:,4)>=1 & ts_id(:,4)<(duration+1) & ismember(ts_id(:,5),trial_sel),:);
                end
                if opt.shuf_trl
                    ts_id_shuf=[];
                    for in_chain_pos=1:numel(cids)
                        nodesel=ts_id(:,3)==in_chain_pos;
                        trl_onset=arrayfun(@(x) trials(x,1),ts_id(nodesel,5));
                        delta_ts=ts_id(nodesel,1)-trl_onset;

                        shufmap=containers.Map(num2cell(trial_sel),num2cell(randsample(trial_sel,numel(trial_sel))));
                        shuf_trl=cell2mat(shufmap.values(num2cell(ts_id(nodesel,5))));
                        shuf_onset=arrayfun(@(x) trials(x,1),shuf_trl);
                        ts_id_shuf=[ts_id_shuf;sortrows([shuf_onset+delta_ts,ts_id(nodesel,2:end)],1)];
                    end
                
                    ts=chain_alt(ts_id_shuf(:,[1 3]));
                else
                    if opt.anti_dir
                        anti_ts_id=ts_id(:,[1 3]);
                        anti_ts_id(:,2)=(numel(cids)+1)-anti_ts_id(:,2);
                        ts=chain_alt(anti_ts_id);
                    else
                        ts=chain_alt(ts_id(:,[1 3]));
                    end
                end

                if ~isempty(ts)
                    outkey="s"+sessid+"c"+cc;
                    out.("d"+duration).(outid).(outkey).ts=ts;
                    out.("d"+duration).(outid).(outkey).meta={cids,chains.tcoms(cc)};
                    if opt.skip_ts_id
                        out.("d"+duration).(outid).(outkey).trials=trials;
                    else
                        out.("d"+duration).(outid).(outkey).ts_id=table(uint32(ts_id(:,1)),uint16(ts_id(:,2)),uint8(ts_id(:,3)),ts_id(:,4),uint16(ts_id(:,5)),'VariableNames',{'TS','CID','POS','Time','Trial'});
                    end

                    % CCG
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
                        out.("d"+duration).(outid).(outkey).ccgs=chainccg;
                    end
                else
                    notfound=[notfound;{sessid},{cc},{cids}];
                end
            end
        end
    end
end
if opt.ccg
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
end

if ~opt.skip_save
    blame=vcs.blame();
    if opt.rev
        save(fullfile('bzdata','chain_rev_tag.mat'),'out','blame');
    elseif opt.shuf_trl
        save(fullfile("bzdata","chain_shuf_tag_"+opt.shuf_idx+".mat"),'out','blame');
    else
        save(fullfile('bzdata','chain_tag.mat'),'out','blame');
    end
end

% stats(out);
end
% 
% function stats(out)
% for dur=reshape(fieldnames(out),1,[])
%     perchaindur=struct();
%     [perchaindur.size,perchaindur.dur]=deal([]);
%     for wv=reshape(fieldnames(out.(dur{1})),1,[])
%         for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
%             perchaindur.size=[perchaindur.size;size(out.(dur{1}).(wv{1}).(lp{1}).ts,2)];
%             perchaindur.dur=[perchaindur.dur;{diff(out.(dur{1}).(wv{1}).(lp{1}).ts(:,[1,end]),1,2)./30}];
%         end
%     end
%     statss.("d"+dur)=perchaindur;
% end
% 
% disp([max(cell2mat(statss.dd3.dur)),max(cell2mat(statss.dd6.dur))])
% 
% end


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


function replay(sschain_trl)
for dd=reshape(fieldnames(sschain_trl),1,[])
    for ww=reshape(fieldnames(sschain_trl.(dd{1})),1,[])
        for cc=reshape(fieldnames(sschain_trl.(dd{1}).(ww{1})),1,[])    
            sessid=str2double(regexp(cc{1},'(?<=s)\d{1,3}(?=c)','match','once'));
            %TODO: move to pre-processing?
            [~,SPKTS,~,~,~,~]=ephys.getSPKID_TS(sessid,'keep_trial',false);
    

            onechain=sschain_trl.(dd{1}).(ww{1}).(cc{1});
            
            dur_pref=str2double(dd{1}(2:end));
            if contains(ww,'s1')
                pref_trl=sschain_trl.(dd{1}).(ww{1}).(cc{1}).trials(:,5)==4 & sschain_trl.(dd{1}).(ww{1}).(cc{1}).trials(:,8)==dur_pref;
            elseif contains(ww,'s2')
                pref_trl=sschain_trl.(dd{1}).(ww{1}).(cc{1}).trials(:,5)==8 & sschain_trl.(dd{1}).(ww{1}).(cc{1}).trials(:,8)==dur_pref;
            else
                keyboard()
            end
            
            sps=30000;

            % [nearest before; nearest after] * [trl_id,dT, samp, delay,performace, wt, prefer]% [nearest before; nearest after] * [trl_id,dT, samp, delay,performace, wt, prefer]
            trl_align=nan(size(onechain.ts,1),14);
            
            for mii=1:size(onechain.ts,1)
                nxt_trl=find(onechain.trials(:,1)>onechain.ts(mii,1),1,"first");
                if nxt_trl==1 % before first
                    trl_align(mii,:)=[-1,-1,-1,-1,-1,-1,-1,nxt_trl,(onechain.trials(nxt_trl,1)-onechain.ts(mii,1))./sps,onechain.trials(nxt_trl,[5 8 9 10]),pref_trl(nxt_trl)];
                elseif isempty(nxt_trl) % after last
                    prev_trl=size(onechain.trials,1);
                    trl_align(mii,:)=[prev_trl,(onechain.ts(mii,1)-onechain.trials(prev_trl,1))./sps,onechain.trials(prev_trl,[5 8 9 10]),pref_trl(prev_trl),-1,-1,-1,-1,-1,-1,-1];
                else % in session
                    prev_trl=nxt_trl-1;
                    trl_align(mii,:)=[prev_trl,(onechain.ts(mii,1)-onechain.trials(prev_trl,1))./sps,onechain.trials(prev_trl,[5 8 9 10]),pref_trl(prev_trl),nxt_trl,(onechain.trials(nxt_trl,1)-onechain.ts(mii,1))./sps,onechain.trials(nxt_trl,[5 8 9 10]),pref_trl(nxt_trl)];
                end
            end
            

            % 1:prefer leading ITI, 2:prefer delay, 3:prefer following ITI, 4:nonprefer leading
            % ITI, 5:nonprefer delay, 6:nonprefer following ITI, 7:much  earlier,
            % 8:much later
            len=size(onechain.ts,2);
            lastTrl=size(onechain.trials,1);
            freqstats=struct();
            pref_trl_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            pref_trl_after_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=(trl_align(:,4)+1);
            
            freqstats.pref_delay_correct=[nnz(pref_trl_delay).*len,(sum(trl_align(pref_trl_delay,4)))];
            freqstats.pref_after_delay_correct=[nnz(pref_trl_after_delay).*len,sum(sum(trl_align(pref_trl_after_delay,[2 9]),2)-trl_align(pref_trl_after_delay,4)-2)];

            lastSps=SPKTS(end);
            freqstats.before_session=[nnz(trl_align(:,8)==1 & trl_align(:,9)>60).*len,(onechain.trials(1,1)./sps-60)];
            freqstats.after_session=[nnz(trl_align(:,1)==lastTrl & trl_align(:,2)>(60+2+trl_align(:,4))).*len,((lastSps-onechain.trials(end,2))./sps-60-1)];

            sschain_trl.(dd{1}).(ww{1}).(cc{1}).trl_align=trl_align;
            sschain_trl.(dd{1}).(ww{1}).(cc{1}).freqstats=freqstats;
        end
    end
end








end
