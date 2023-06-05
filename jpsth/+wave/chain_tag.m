%
% load(fullfile('bzdata','chains_mix.mat'),'chains_uf');
% chains=chains_uf;
% clear chains_uf;
classdef chain_tag < handle
    methods(Static)
        function [out,notfound]=tag(chains,opt)
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

                                ts=wave.chain_tag.chain_alt(ts_id_shuf(:,[1 3]));
                            else
                                if opt.anti_dir
                                    anti_ts_id=ts_id(:,[1 3]);
                                    anti_ts_id(:,2)=(numel(cids)+1)-anti_ts_id(:,2);
                                    ts=wave.chain_tag.chain_alt(anti_ts_id);
                                else
                                    ts=wave.chain_tag.chain_alt(ts_id(:,[1 3]));
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


        function sschain_trl=replay(sschain_trl,ttl)
            arguments
                sschain_trl
                ttl (1,:) char = 'chains'
            end
            for dd=reshape(fieldnames(sschain_trl),1,[])
                for ww=reshape(fieldnames(sschain_trl.(dd{1})),1,[])
                    for cc=reshape(fieldnames(sschain_trl.(dd{1}).(ww{1})),1,[])
                        sessid=str2double(regexp(cc{1},'(?<=s)\d{1,3}(?=c)','match','once'));
                        %TODO: move to pre-processing?
                        [~,SPKTS,~,~,~,~]=ephys.getSPKID_TS(sessid,'keep_trial',false);

                        onechain=sschain_trl.(dd{1}).(ww{1}).(cc{1});

                        dur_pref=str2double(dd{1}(2:end));
                        if contains(ww,'s1')
                            samp_pref=4;
                        elseif contains(ww,'s2')
                            samp_pref=8;
                        else
                            disp('no wave id')
                            keyboard()
                        end

                        pref_trl=sschain_trl.(dd{1}).(ww{1}).(cc{1}).trials(:,5)==samp_pref & sschain_trl.(dd{1}).(ww{1}).(cc{1}).trials(:,8)==dur_pref;

                        sps=30000;

                        % [nearest before; nearest after] * [trl_id,dT, samp, delay,wt,correct,prefer]% [nearest before; nearest after] * [trl_id,dT, samp, delay,performace, wt, prefer]
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

                        len=size(onechain.ts,2);
                        lastTrl=size(onechain.trials,1);
                        freqstats=struct();

                        % delay correct
                        pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                        pref_delay_trls=all(onechain.trials(:,9:10)==1,2) & onechain.trials(:,5)==samp_pref & onechain.trials(:,8)==dur_pref;
                        freqstats.pref_delay_correct=[nnz(pref_delay).*len,sum(onechain.trials(pref_delay_trls,8))];

                        % delay error
                        pref_delay_err=trl_align(:,6)==0 & trl_align(:,7)==1 & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                        pref_delay_err_trls=onechain.trials(:,10)==0 & onechain.trials(:,5)==samp_pref & onechain.trials(:,8)==dur_pref;
                        freqstats.pref_delay_error=[nnz(pref_delay_err).*len,sum(onechain.trials(pref_delay_err_trls,8))];

                        % delay non-prefered
                        nonpref_delay=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp_pref & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                        nonpref_delay_trls=all(onechain.trials(:,9:10)==1,2) & onechain.trials(:,5)~=samp_pref;
                        freqstats.nonpref_delay_correct=[nnz(nonpref_delay).*len,(sum(onechain.trials(nonpref_delay_trls,8)))];

                        % 1/samp delay 1/test 2/rwd?
                        % decision/test correct
                        pref_test=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=(trl_align(:,4)+1) & trl_align(:,2)<(trl_align(:,4)+2);
                        freqstats.pref_test=[nnz(pref_test).*len,nnz(pref_delay_trls)]; % 1 sec per trl

                        %supply for last trial
                        onechain.trials(end+1,:)=onechain.trials(end,2)+14*sps;
                        
                        % succeed ITI pref correct
                        pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                            & trl_align(:,2)>=(trl_align(:,4)+4)...  % not in delay or test
                            & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
                        freqstats.pref_succeed_ITI=[nnz(pref_succeed_iti).*len,...
                            sum((onechain.trials(find(pref_delay_trls)+1,1)-onechain.trials(pref_delay_trls,2))./sps-3)]; %rwd + test

                        % succeed ITI nonpref
                        nonpref_succeed_iti=(all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp_pref) ... % WT, nonpref
                            & trl_align(:,2)>=(trl_align(:,4)+4)...  % not in delay or test
                            & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14)); % 1s samp, 1s test, 2s rwd
                        freqstats.nonpref_succeed_ITI=[nnz(nonpref_succeed_iti).*len,...
                            sum((onechain.trials(find(nonpref_delay_trls)+1,1)-onechain.trials(nonpref_delay_trls,2))./sps-3)];

                        onechain.trials(end,:)=[];

                        % TODO: precede preferred, non preferred
%                       
%                         % precede ITI pref correct
                        onechain.trials=[repmat(onechain.trials(1,1)-14*sps,1,14);onechain.trials];
                        pref_precede_iti=all(trl_align(:,12:14)==1,2) & (trl_align(:,2)>=(trl_align(:,4)+4) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
                        freqstats.pref_precede_ITI=[nnz(pref_precede_iti).*len,sum((onechain.trials(pref_delay_trls,1)-onechain.trials(find(pref_delay_trls)-1,2))./sps-3)]; %rwd + test
% 
%                         %  precede ITI nonpref
%                         nonpref_late_iti=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp_pref & trl_align(:,2)>=(trl_align(:,4)+4); % 1s samp, 1s test, 2s rwd
%                         freqstats.nonpref_late_ITI=[nnz(nonpref_late_iti).*len,sum((onechain.trials(find(nonpref_delay_trls)+1,1)-onechain.trials(nonpref_delay_trls,2))./sps-3)];



                        % long before and after 
                        lastSps=SPKTS(end);
                        freqstats.before_session=[nnz(trl_align(:,8)==1 & trl_align(:,9)>60).*len,onechain.trials(1,1)./sps-60];
                        freqstats.after_session=[nnz(trl_align(:,1)==lastTrl & trl_align(:,2)>(60+2+trl_align(:,4))).*len,(lastSps-onechain.trials(end,2))./sps-60-1];

                        sschain_trl.(dd{1}).(ww{1}).(cc{1}).trl_align=trl_align;
                        sschain_trl.(dd{1}).(ww{1}).(cc{1}).freqstats=freqstats;
                    end
                end
            end

            stats=[];
            for dd=reshape(fieldnames(sschain_trl),1,[])
                for ww=reshape(fieldnames(sschain_trl.(dd{1})),1,[])
                    for cc=reshape(fieldnames(sschain_trl.(dd{1}).(ww{1})),1,[])
                        stats=[stats,cellfun(@(x) x(1)./x(2),struct2cell(sschain_trl.(dd{1}).(ww{1}).(cc{1}).freqstats))];
                    end
                end
            end

            figure();
            boxplot(stats.','Colors','k','Symbol','c.')
            ylim([0,1.5])
            set(gca(),'XTick',1:7,'XTickLabel',fieldnames(sschain_trl.(dd{1}).(ww{1}).(cc{1}).freqstats),'YScale','linear')
            for jj=2:7
                pp=ranksum(stats(1,:),stats(jj,:));
                text(jj,1.5,sprintf('%.3f',pp),'VerticalAlignment','top','HorizontalAlignment','center')
            end
            title(ttl)

            figure();
            hold on
            for ii=1:7
                swarmchart(repmat(ii,size(stats,2),1),stats(ii,:),16,'.')
            end
            ylim([-0.1,1.5])
            set(gca(),'XTick',1:7,'XTickLabel',fieldnames(sschain_trl.(dd{1}).(ww{1}).(cc{1}).freqstats),'YScale','linear','TickLabelInterpreter','none')

        end
    end
end