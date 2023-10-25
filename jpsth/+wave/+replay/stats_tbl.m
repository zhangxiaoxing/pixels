function [chains,loops]=stats_tbl(sschain_trl,ssloop_trl,trials_dict,opt)
arguments
    sschain_trl = []
    ssloop_trl = []
    trials_dict = []
    opt.var_len (1,1) logical = false
    opt.cross_only (1,1) logical = false
    opt.within_only (1,1) logical = false
    opt.nonmem (1,1) logical = false
    opt.skip_save (1,1) logical = true
    % opt.nonmem_ring (1,1) logical = false
    opt.shuf (1,1) logical = false
    opt.shufidx = 1
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
assert(~(opt.cross_only && opt.within_only),"conflict selection")

if isempty(trials_dict)
    switch opt.criteria
        case 'WT'
            load(fullfile('binary','trials_dict.mat'),'trials_dict');
        case 'Learning'
            trials_dict=behav.get_trials_dict('skip_save',true,'criteria','Learning');
        otherwise
            keyboard()
    end
end

if isempty(sschain_trl)
    if opt.shuf
        % TODO: LN WIP
        switch opt.criteria
            case 'WT'
                fstr=load(fullfile('binary','shufs',sprintf('chain_tag_shuf%d.mat',opt.shufidx)),'out');
            case 'Learning'
                fstr=load(fullfile('binary','shufs',sprintf('LN_chain_tag_shuf%d.mat',opt.shufidx)),'out');
            otherwise
                keyboard();
        end
    elseif opt.nonmem
        error("Learning not ready")
        fstr=load(fullfile('binary','chain_tag_nonmem_all_trl.mat'),'out');
        opt.var_len=false;
        [chain_replay,chain_sums,chain_raw]=stats_one(fstr.out,trials_dict,opt);
        blame=vcs.blame();
        save(fullfile('binary','motif_replay_chain_nonmem.mat'),'chain_raw','chain_sums','chain_replay','-v7.3');
    else
        switch opt.criteria
            case 'WT'
                fstr=load(fullfile('binary','chain_tag_all_trl.mat'),'out');
            case 'Learning'
                fstr=load(fullfile('binary','LN_chain_tag_all_trl.mat'),'out');
            otherwise
                keyboard();
        end
    end
    sschain_trl=fstr.out;
    clear fstr
end

if isempty(ssloop_trl)
    if opt.shuf
        switch opt.criteria
            case 'WT'
                load(fullfile('binary','shufs',sprintf('ring_tag_shuf%d.mat',opt.shufidx)),'ssloop_trl');
            case 'Learning'
                load(fullfile('binary','shufs',sprintf('LN_ring_tag_shuf%d.mat',opt.shufidx)),'ssloop_trl');
            otherwise
                keyboard();
        end
    elseif opt.nonmem
        error("Not tested yet")
        switch opt.criteria
            case 'WT'
                load(fullfile('binary','rings_tag_trl.mat'),'ssloop_trl')
            case 'Learning'
                load(fullfile('binary','LN_rings_tag_trl.mat'),'ssloop_trl')
            otherwise
                keyboard();
        end
        opt.var_len=true;
        [ring_replay,loops_sums,loops_raw]=stats_one(ssloop_trl,trials_dict,opt);
        blame=vcs.blame();
        save(fullfile('binary','motif_replay_ring_nonmem.mat'),'loops_raw','loops_sums','ring_replay','blame','opt','-v7.3');
    else
        switch opt.criteria
            case 'WT'
                load(fullfile('binary','rings_tag_trl.mat'),'ssloop_trl')
            case 'Learning'
                load(fullfile('binary','LN_rings_tag_trl.mat'),'ssloop_trl')
            otherwise
                keyboard();
        end
    end
end


opt.var_len=false;
[chain_replay,chain_sums,chain_raw]=stats_one(sschain_trl,trials_dict,opt);
opt.var_len=true;
[ring_replay,loops_sums,loops_raw]=stats_one(ssloop_trl,trials_dict,opt);

chains=cell2struct({chain_replay;chain_sums;chain_raw},{'replay','sums','raw'});
loops=cell2struct({ring_replay;loops_sums;loops_raw},{'replay','sums','raw'});

if ~opt.skip_save
    blame=vcs.blame();
    if opt.shuf
        switch opt.criteria
            case 'WT'
                save(fullfile("binary","motif_replay_shuf"+opt.shufidx+".mat"),'chain_replay','ring_replay','chain_sums','loops_sums','chain_raw','loops_raw','blame','opt');
            case 'Learning'
                save(fullfile("binary","shufs","LN_motif_replay_shuf"+opt.shufidx+".mat"),'chain_replay','ring_replay','chain_sums','loops_sums','chain_raw','loops_raw','blame','opt');
            otherwise
                keyboard();
        end

    else
        switch opt.criteria
            case 'WT'
                save(fullfile("binary","motif_replay.mat"),'chain_replay','ring_replay','chain_sums','loops_sums','chain_raw','loops_raw','blame','opt','-v7.3');
            case 'Learning'
                save(fullfile("binary","LN_motif_replay.mat"),'chain_replay','ring_replay','chain_sums','loops_sums','chain_raw','loops_raw','blame','opt','-v7.3');
            otherwise
                keyboard();
        end
    end
end
end

function [motif_replay,sum_stats,raw]=stats_one(motif_replay,trials_dict,opt)
sps=30000;
nmsel=(motif_replay.wave=="none" | contains(motif_replay.wave,'nm'));
if opt.nonmem
    motif_replay=motif_replay(nmsel,:);
else
    motif_replay=motif_replay(~nmsel,:);
end
    
stat_cell=cell(size(motif_replay,1),2);
for tidx=1:size(motif_replay,1)
    % onechain=motif_replay.(dd{1}).(ww{1}).(cc{1});
    if (opt.cross_only && ~motif_replay.meta{tidx,3}) ...
            || (opt.within_only && motif_replay.meta{tidx,3})
        continue
    end

    trials=cell2mat(trials_dict(motif_replay.session(tidx)));
    if strcmp(opt.criteria,'Learning')
        trials=behav.procPerf(trials,"criteria","Learning");
        trials(:,10)=trials(:,9);
    end
    session_tick=wave.replay.sessid2length(motif_replay.session(tidx),'criteria',opt.criteria);
    if ~(strcmp(motif_replay.wave(tidx),"none") || contains(motif_replay.wave(tidx),'nm'))
        dur_pref=motif_replay.delay(tidx);
        if contains(motif_replay.wave(tidx),"s1")
            samp_pref=4;
        elseif contains(motif_replay.wave(tidx),"s2")
            samp_pref=8;
        else
            disp('no wave id')
            keyboard()
        end
        pref_trl=trials(:,5)==samp_pref & trials(:,8)==dur_pref;
        nnonmem=false;
    else %non mem
        pref_trl=ismember(trials(:,5),[4 8]) & ismember(trials(:,8),[3 6]);
        nnonmem=true;
    end

    % [nearest before; nearest after] * [trl_id,dT, samp, delay,wt,correct,prefer]% [nearest before; nearest after] * [trl_id,dT, samp, delay,performace, wt, prefer]
    if opt.var_len && ~iscell(motif_replay.ts{tidx})
        motif_replay.ts{tidx}=motif_replay.ts(tidx);
    end

    trl_align=nan(size(motif_replay.ts{tidx},1),14);

    for mii=1:size(motif_replay.ts{tidx},1)
        if opt.var_len
            one_onset=motif_replay.ts{tidx}{mii}(1);
        else
            one_onset=motif_replay.ts{tidx}(mii,1);
        end
        nxt_trl=find(trials(:,1)>one_onset,1,"first");
        if nxt_trl==1 % before first
            trl_align(mii,:)=[-1,-1,-1,-1,-1,-1,-1,nxt_trl,(trials(nxt_trl,1)-one_onset)./sps,trials(nxt_trl,[5 8 9 10]),pref_trl(nxt_trl)];
        elseif isempty(nxt_trl) % after last
            prev_trl=size(trials,1);
            trl_align(mii,:)=[prev_trl,(one_onset-trials(prev_trl,1))./sps,trials(prev_trl,[5 8 9 10]),pref_trl(prev_trl),-1,-1,-1,-1,-1,-1,-1];
        else % in session
            prev_trl=nxt_trl-1;
            trl_align(mii,:)=[prev_trl,(one_onset-trials(prev_trl,1))./sps,trials(prev_trl,[5 8 9 10]),pref_trl(prev_trl),nxt_trl,(trials(nxt_trl,1)-one_onset)./sps,trials(nxt_trl,[5 8 9 10]),pref_trl(nxt_trl)];
        end
    end
    lastTrl=size(trials,1);
    freqstats=struct();

    if opt.var_len
        len=cellfun(@(x) numel(x),motif_replay.ts{tidx});
    else
        len=size(motif_replay.ts{tidx},2);
    end

    % delay correct
    pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
    if nnonmem
        pref_delay_trls=all(trials(:,9:10)==1,2) & ismember(trials(:,5),[4 8]) & ismember(trials(:,8),[3 6]);
    else
        pref_delay_trls=all(trials(:,9:10)==1,2) & trials(:,5)==samp_pref & trials(:,8)==dur_pref;
    end
    freqstats.pref_delay_correct=[sum(pref_delay.*len),sum(trials(pref_delay_trls,8))];
    
    if ~nnonmem
        % delay error
        pref_delay_err=trl_align(:,6)==0 & trl_align(:,7)==1 & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
        pref_delay_err_trls=trials(:,10)==0 & trials(:,5)==samp_pref & trials(:,8)==dur_pref;
        freqstats.pref_delay_error=[sum(pref_delay_err.*len),sum(trials(pref_delay_err_trls,8))];

        % delay non-prefered
        nonpref_delay=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp_pref & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
        nonpref_delay_trls=all(trials(:,9:10)==1,2) & trials(:,5)~=samp_pref;
        freqstats.nonpref_delay_correct=[sum(nonpref_delay.*len),(sum(trials(nonpref_delay_trls,8)))];

        % 1/samp delay 1/test 2/rwd?
        % decision/test correct
        pref_test=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=(trl_align(:,4)+1) & trl_align(:,2)<(trl_align(:,4)+2);
        freqstats.pref_test=[sum(pref_test.*len),nnz(pref_delay_trls)]; % 1 sec per trl
    end
    %supply for last trial
    trials(end+1,:)=min(session_tick-3,trials(end,2)+14*sps);

    % succeed ITI pref correct
    pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
        & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test
        & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
    freqstats.pref_succeed_ITI=[sum(pref_succeed_iti.*len),...
        sum((trials(find(pref_delay_trls)+1,1)-trials(pref_delay_trls,2))./sps-3)]; %rwd + test

    if ~nnonmem
        % succeed ITI pref error
        pref_succeed_iti_err=trl_align(:,6)==0 & trl_align(:,7)==1 ...
            & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test
            & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
        freqstats.pref_succeed_ITI_err=[sum(pref_succeed_iti_err.*len),...
            sum((trials(find(pref_delay_err_trls)+1,1)-trials(pref_delay_err_trls,2))./sps-3)]; %rwd + test

        % succeed ITI nonpref
        nonpref_succeed_iti=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp_pref ... % WT, nonpref
            & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test
            & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14)); % 1s samp, 1s test, 2s rwd
        freqstats.nonpref_succeed_ITI=[sum(nonpref_succeed_iti.*len),...
            sum((trials(find(nonpref_delay_trls)+1,1)-trials(nonpref_delay_trls,2))./sps-3)];
    end
    trials(end,:)=[];

    % precede preferred, non preferred

    % precede ITI pref correct
    trials=[repmat(trials(1,1)-14*sps,1,10);trials];
    pref_precede_iti=all(trl_align(:,12:14)==1,2)...
        & (trl_align(:,2)>=(trl_align(:,4)+5) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
    freqstats.pref_precede_ITI=[sum(pref_precede_iti.*len),sum((trials(find(pref_delay_trls)+1,1)-trials(pref_delay_trls,2))./sps-3)]; %rwd + test
    if ~nnonmem
        % precede ITI pref error
        pref_succeed_iti_err=trl_align(:,13)==0 & trl_align(:,14)==1 ...
            & (trl_align(:,2)>=(trl_align(:,4)+5) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
        freqstats.pref_precede_ITI_err=[sum(pref_succeed_iti_err.*len),sum((trials(find(pref_delay_err_trls)+1,1)-trials(pref_delay_err_trls,2))./sps-3)]; %rwd + test

        %  precede ITI nonpref
        nonpref_precede_iti=all(trl_align(:,12:13)==1,2) & trl_align(:,10)~=samp_pref...
            & (trl_align(:,2)>=(trl_align(:,4)+5) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
        freqstats.nonpref_precede_ITI=[sum(nonpref_precede_iti.*len),sum((trials(find(nonpref_delay_trls)+1,1)-trials(nonpref_delay_trls,2))./sps-3)]; %rwd + test
    end
    trials(1,:)=[];

    % long before and after
    sessid=motif_replay.session(tidx);
    rec_dur=wave.replay.sessid2length(sessid);
    freqstats.before_session=[sum((trl_align(:,8)==1 & trl_align(:,9)>60).*len),trials(1,1)./sps-60];
    freqstats.after_session=[sum((trl_align(:,1)==lastTrl & trl_align(:,2)>(60+2+trl_align(:,4))).*len),(rec_dur-trials(end,2))./sps-60-1];

    % time length criteria since 23-Jun-28
    if freqstats.before_session(2)<=30
        freqstats.before_session=nan(1,2);
    end

    if freqstats.after_session(2)<=30
        freqstats.after_session=nan(1,2);
    end

    if ~nnonmem
        % nonpreferred error, last-minute addition
        % prev_trl,onset,trials(prev_trl,[5 8 9 10]),pref_trl(prev_trl)
        np_delay_err=trl_align(:,6)==0 & trl_align(:,3)==setdiff([4,8],samp_pref) & trl_align(:,4)==dur_pref & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
        np_delay_err_trls=trials(:,10)==0 & trials(:,5)==setdiff([4,8],samp_pref) & trials(:,8)==dur_pref;
        freqstats.np_delay_error=[sum(np_delay_err.*len),sum(trials(np_delay_err_trls,8))];
    end

    stat_cell(tidx,1)={trl_align};
    stat_cell(tidx,2)={freqstats};
end
motif_replay.trl_align=stat_cell(:,1);
motif_replay.freqstats=stat_cell(:,2);

sum_stats=[];
raw=cell2struct({[];[];cell(0);cell(0)},{'count','time','condition','tag'});

motiftype=motif_replay.Properties.VariableNames{4};

for tidx=1:size(motif_replay,1)
    if (opt.cross_only && ~motif_replay.meta{tidx,3}) ...
            || (opt.within_only && motif_replay.meta{tidx,3})
        continue
    end

    sum_stats=[sum_stats,cellfun(@(x) x(1)./x(2),struct2cell(motif_replay.freqstats{tidx}))];
    raw.count=[raw.count,cellfun(@(x) x(1),struct2cell(motif_replay.freqstats{tidx}))];
    raw.time=[raw.time,cellfun(@(x) x(2),struct2cell(motif_replay.freqstats{tidx}))];
    sess_str=regexp(motif_replay.(motiftype)(tidx),'^s\d{1,3}(?=(c|r))','match','once');
    raw.condition{end+1}=sess_str+replace(motif_replay.wave(tidx),"s","o");
    raw.tag{end+1}=regexp(motif_replay.(motiftype)(tidx),'(c|r)\d*$','match','once');
end

% fn=replace(motiftype,"_id","")+"_replay.mat";
% TODO: streamline


end

