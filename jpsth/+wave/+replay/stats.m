function [motif_replay,sum_stats,raw]=stats(motif_replay,opt)
arguments
    motif_replay
    opt.var_len (1,1) logical = false
    opt.cross_only (1,1) logical = false
    opt.within_only (1,1) logical = false
    % opt.nonmem_ring (1,1) logical = false
end
assert(~(opt.cross_only && opt.within_only),"conflict selection")


sps=30000;

for dd=reshape(fieldnames(motif_replay),1,[])
    for ww=reshape(fieldnames(motif_replay.(dd{1})),1,[])
        for cc=reshape(fieldnames(motif_replay.(dd{1}).(ww{1})),1,[])
            onechain=motif_replay.(dd{1}).(ww{1}).(cc{1});
            if (opt.cross_only && ~onechain.meta{3}) ...
                    || (opt.within_only && onechain.meta{3})
                continue
            end

            if ~(strcmp(ww{1},'none') || contains(ww,'nm'))
                dur_pref=str2double(dd{1}(2:end));
                if contains(ww,'s1')
                    samp_pref=4;
                elseif contains(ww,'s2')
                    samp_pref=8;
                else
                    disp('no wave id')
                    keyboard()
                end
                pref_trl=onechain.trials(:,5)==samp_pref & onechain.trials(:,8)==dur_pref;
                nnonmem=false;
            else %non mem
                pref_trl=ismember(onechain.trials(:,5),[4 8]) & ismember(onechain.trials(:,8),[3 6]);
                nnonmem=true;
            end

            % [nearest before; nearest after] * [trl_id,dT, samp, delay,wt,correct,prefer]% [nearest before; nearest after] * [trl_id,dT, samp, delay,performace, wt, prefer]
            trl_align=nan(size(onechain.ts,1),14);

            for mii=1:size(onechain.ts,1)
                if opt.var_len
                    one_onset=onechain.ts{mii}(1);
                else
                    one_onset=onechain.ts(mii,1);
                end
                nxt_trl=find(onechain.trials(:,1)>one_onset,1,"first");
                if nxt_trl==1 % before first
                    trl_align(mii,:)=[-1,-1,-1,-1,-1,-1,-1,nxt_trl,(onechain.trials(nxt_trl,1)-one_onset)./sps,onechain.trials(nxt_trl,[5 8 9 10]),pref_trl(nxt_trl)];
                elseif isempty(nxt_trl) % after last
                    prev_trl=size(onechain.trials,1);
                    trl_align(mii,:)=[prev_trl,(one_onset-onechain.trials(prev_trl,1))./sps,onechain.trials(prev_trl,[5 8 9 10]),pref_trl(prev_trl),-1,-1,-1,-1,-1,-1,-1];
                else % in session
                    prev_trl=nxt_trl-1;
                    trl_align(mii,:)=[prev_trl,(one_onset-onechain.trials(prev_trl,1))./sps,onechain.trials(prev_trl,[5 8 9 10]),pref_trl(prev_trl),nxt_trl,(onechain.trials(nxt_trl,1)-one_onset)./sps,onechain.trials(nxt_trl,[5 8 9 10]),pref_trl(nxt_trl)];
                end
            end
            lastTrl=size(onechain.trials,1);
            freqstats=struct();

            if opt.var_len
                len=cellfun(@(x) numel(x),onechain.ts);
            else
                len=size(onechain.ts,2);
            end

            % delay correct
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            if nnonmem
                pref_delay_trls=all(onechain.trials(:,9:10)==1,2) & ismember(onechain.trials(:,5),[4 8]) & ismember(onechain.trials(:,8),[3 6]);
            else
                pref_delay_trls=all(onechain.trials(:,9:10)==1,2) & onechain.trials(:,5)==samp_pref & onechain.trials(:,8)==dur_pref;
            end
            freqstats.pref_delay_correct=[sum(pref_delay.*len),sum(onechain.trials(pref_delay_trls,8))];

            % delay error
            if ~nnonmem
                pref_delay_err=trl_align(:,6)==0 & trl_align(:,7)==1 & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                pref_delay_err_trls=onechain.trials(:,10)==0 & onechain.trials(:,5)==samp_pref & onechain.trials(:,8)==dur_pref;
                freqstats.pref_delay_error=[sum(pref_delay_err.*len),sum(onechain.trials(pref_delay_err_trls,8))];

                % delay non-prefered
                nonpref_delay=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp_pref & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                nonpref_delay_trls=all(onechain.trials(:,9:10)==1,2) & onechain.trials(:,5)~=samp_pref;
                freqstats.nonpref_delay_correct=[sum(nonpref_delay.*len),(sum(onechain.trials(nonpref_delay_trls,8)))];

                % 1/samp delay 1/test 2/rwd?
                % decision/test correct
                pref_test=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=(trl_align(:,4)+1) & trl_align(:,2)<(trl_align(:,4)+2);
                freqstats.pref_test=[sum(pref_test.*len),nnz(pref_delay_trls)]; % 1 sec per trl
            end
            %supply for last trial
            onechain.trials(end+1,:)=onechain.trials(end,2)+14*sps;

            % succeed ITI pref correct
            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            freqstats.pref_succeed_ITI=[sum(pref_succeed_iti.*len),...
                sum((onechain.trials(find(pref_delay_trls)+1,1)-onechain.trials(pref_delay_trls,2))./sps-3)]; %rwd + test

            if ~nnonmem
                % succeed ITI pref error
                pref_succeed_iti_err=trl_align(:,6)==0 & trl_align(:,7)==1 ...
                    & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
                freqstats.pref_succeed_ITI_err=[sum(pref_succeed_iti_err.*len),...
                    sum((onechain.trials(find(pref_delay_err_trls)+1,1)-onechain.trials(pref_delay_err_trls,2))./sps-3)]; %rwd + test

                % succeed ITI nonpref
                nonpref_succeed_iti=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp_pref ... % WT, nonpref
                    & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14)); % 1s samp, 1s test, 2s rwd
                freqstats.nonpref_succeed_ITI=[sum(nonpref_succeed_iti.*len),...
                    sum((onechain.trials(find(nonpref_delay_trls)+1,1)-onechain.trials(nonpref_delay_trls,2))./sps-3)];
            end
            onechain.trials(end,:)=[];

            % precede preferred, non preferred

            % precede ITI pref correct
            onechain.trials=[repmat(onechain.trials(1,1)-14*sps,1,10);onechain.trials];
            pref_precede_iti=all(trl_align(:,12:14)==1,2)...
                & (trl_align(:,2)>=(trl_align(:,4)+5) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
            freqstats.pref_precede_ITI=[sum(pref_precede_iti.*len),sum((onechain.trials(find(pref_delay_trls)+1,1)-onechain.trials(pref_delay_trls,2))./sps-3)]; %rwd + test
            if ~nnonmem
                % precede ITI pref error
                pref_succeed_iti_err=trl_align(:,13)==0 & trl_align(:,14)==1 ...
                    & (trl_align(:,2)>=(trl_align(:,4)+5) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
                freqstats.pref_precede_ITI_err=[sum(pref_succeed_iti_err.*len),sum((onechain.trials(find(pref_delay_err_trls)+1,1)-onechain.trials(pref_delay_err_trls,2))./sps-3)]; %rwd + test

                %  precede ITI nonpref
                nonpref_precede_iti=all(trl_align(:,12:13)==1,2) & trl_align(:,10)~=samp_pref...
                    & (trl_align(:,2)>=(trl_align(:,4)+5) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
                freqstats.nonpref_precede_ITI=[sum(nonpref_precede_iti.*len),sum((onechain.trials(find(nonpref_delay_trls)+1,1)-onechain.trials(nonpref_delay_trls,2))./sps-3)]; %rwd + test
            end
            onechain.trials(1,:)=[];
            %
            % if nnz(nonpref_precede_iti)==0 && nnz(nonpref_succeed_iti)>0
            %     keyboard()
            % end

            % long before and after
            sessid=str2double(regexp(cc{1},'(?<=s)\d{1,3}(?=(c|r))','match','once'));
            rec_dur=wave.replay.sessid2length(sessid);
            freqstats.before_session=[sum((trl_align(:,8)==1 & trl_align(:,9)>60).*len),onechain.trials(1,1)./sps-60];
            freqstats.after_session=[sum((trl_align(:,1)==lastTrl & trl_align(:,2)>(60+2+trl_align(:,4))).*len),(rec_dur-onechain.trials(end,2))./sps-60-1];

            % time length criteria since 23-Jun-28
            if freqstats.before_session(2)<=30
                freqstats.before_session=nan(1,2);
            end

            if freqstats.after_session(2)<=30
                freqstats.after_session=nan(1,2);
            end

            motif_replay.(dd{1}).(ww{1}).(cc{1}).trl_align=trl_align;
            motif_replay.(dd{1}).(ww{1}).(cc{1}).freqstats=freqstats;
        end
    end
end

sum_stats=[];
raw=cell2struct({[];[];cell(0);cell(0)},{'count','time','condition','tag'});
for dd=reshape(fieldnames(motif_replay),1,[])
    for ww=reshape(fieldnames(motif_replay.(dd{1})),1,[])
        for cc=reshape(fieldnames(motif_replay.(dd{1}).(ww{1})),1,[])

            if (opt.cross_only && ~motif_replay.(dd{1}).(ww{1}).(cc{1}).meta{3}) ...
                    || (opt.within_only && motif_replay.(dd{1}).(ww{1}).(cc{1}).meta{3})
                continue
            end

            sum_stats=[sum_stats,cellfun(@(x) x(1)./x(2),struct2cell(motif_replay.(dd{1}).(ww{1}).(cc{1}).freqstats))];
            raw.count=[raw.count,cellfun(@(x) x(1),struct2cell(motif_replay.(dd{1}).(ww{1}).(cc{1}).freqstats))];
            raw.time=[raw.time,cellfun(@(x) x(2),struct2cell(motif_replay.(dd{1}).(ww{1}).(cc{1}).freqstats))];
            sess_str=regexp(cc{1},'^s\d{1,3}(?=(c|r))','match','once');
            raw.condition{end+1}=[sess_str,replace(ww{1},'olf_s','o'),dd{1}];
            raw.tag{end+1}=regexp(cc{1},'(c|r)\d*$','match','once');
        end
    end
end
end