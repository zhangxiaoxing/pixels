% tag spikes in neuron with loop or chain activity

skipfile=false;

if true%~exist('inited','var') || ~inited  % denovo data generation
    %% single spike chain
    sschain=load('chain_tag.mat','out');
    keys=[struct2cell(structfun(@(x) fieldnames(x), sschain.out.d6, 'UniformOutput', false));...
        struct2cell(structfun(@(x) fieldnames(x), sschain.out.d3, 'UniformOutput', false))];
    keys=vertcat(keys{:});
    ssc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));

    %% single spike loop
    load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats'); % bz.rings.rings_time_constant
    pstats=rmfield(pstats,"nonmem");
    ssl_sess=unique(str2double(regexp(fieldnames(pstats.congru),'(?<=s)\d{1,3}(?=r)','match','once')));

    usess=intersect(ssc_sess,ssl_sess);

    % single spk chn:1, burst spk chn:2, single spk loop:4, burst spk loop:8

    %% per-session entrance
    per_trial_motif_cid=cell(0);
    per_trial_motif_freq=[];
    
    for sessid=reshape(usess,1,[])
        disp(sessid)
        [~,~,trials,~,~,~]=ephys.getSPKID_TS(sessid,'keep_trial',false);

        wtsel=trials(:,9)>0 & trials(:,10)>0 & ismember(trials(:,5),[4 8]) &ismember(trials(:,8),[3 6]);
        per_trial_motif_cid=[per_trial_motif_cid;cell(nnz(wtsel),2)];

        per_trial_motif_freq=[per_trial_motif_freq;...
            repmat(sessid,nnz(wtsel),1),find(wtsel),trials(wtsel,5),trials(wtsel,8),...
            zeros(nnz(wtsel),4)];

        %% single spike chain
        for dur=["d6","d3"]
            dd=str2double(replace(dur,"d",""));
            for wid=reshape(fieldnames(sschain.out.(dur)),1,[])
                for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[])
                    if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                        continue
                    end
                    onechain=sschain.out.(dur).(wid{1}).(cc{1});
                    % ts_id: ts, cid, pos, trial_time, trial

                    % per-trial frequency.
                    tsel=[];
                    if dd==3
                        switch wid{1}
                            case 's1d3'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==3;
                            case 's2d3'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==3;
                            case 'olf_s1'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==3;
                            case 'olf_s2'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==3;
                            case 'dur_d3'
                                tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,3),[4 8]) & per_trial_motif_freq(:,4)==3;
                        end
                    else
                        switch wid{1}
                            case 's1d6'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==6;
                            case 's2d6'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==6;
                            case 'olf_s1'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==6;
                            case 'olf_s2'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==6;
                            case 'dur_d6'
                                tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,3),[4 8]) & per_trial_motif_freq(:,4)==6;
                        end
                    end
                    if ~isempty(tsel)
                        for ttt=reshape(find(tsel),1,[])
                            per_trial_motif_cid{ttt,1}=[per_trial_motif_cid{ttt,1},onechain.meta(1)];
                        end

                    end
                end
            end
        end
        %% single spike loop
        for cc=reshape(fieldnames(pstats.congru),1,[])
            if ~startsWith(cc{1},['s',num2str(sessid),'r'])
                continue
            end
            onechain=pstats.congru.(cc{1});
            % .ts_id(:,6) => loop tag
            % per-trial frequency.
            [pref3,pref6]=bz.rings.preferred_trials_rings(onechain.rstats{4},trials);
            tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,2),[pref3;pref6]);
            % per wave
            if strcmp(onechain.rstats{5},'olf')

            elseif strcmp(onechain.rstats{5},'dur')
            
            end
                for ttt=reshape(find(tsel),1,[])
                    per_trial_motif_cid{ttt,2}=[per_trial_motif_cid{ttt,2},onechain.rstats(3)];
                end
                %check cid in largest network
        end
    end
    if ~skipfile
        blame=vcs.blame();
        save(fullfile("bzdata","SingleSpikeChainedLoopCid.mat"),'per_trial_motif_cid','per_trial_motif_freq','blame')
    end
else
    load(fullfile("bzdata","SingleSpikeChainedLoopCid.mat"),'per_trial_motif_cid','per_trial_motif_freq','blame')
end

processed=cell(0);
disconnected=cell(0);

for tt=1:size(per_trial_motif_cid,1)
    if isempty(per_trial_motif_cid{tt,1}) && isempty(per_trial_motif_cid{tt,2})
        continue
    end

    if numel(per_trial_motif_cid{tt,1})+numel(per_trial_motif_cid{tt,2})==1
        % single chain or loop
        % Logically should be dealt with, but practically no instance
        % found.
        keyboard()
    end

    % build graph network
    edges=categorical(unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',[per_trial_motif_cid{tt,:}],'UniformOutput',false).'),'rows'));
    gh=graph(edges(:,1),edges(:,2));
    conncomp=gh.conncomp();
    % check module in network
    if any(conncomp~=1)
        comps=unique(conncomp);
        counter=zeros(numel(comps),1);
        gnodes=cellfun(@(x) str2double(x),gh.Nodes.Name);
        for mm=[per_trial_motif_cid{tt,:}]
            [~,nidx]=ismember(mm{1},gnodes);
            compidx=unique(conncomp(nidx));
            if numel(compidx)==1
                counter(compidx)=counter(compidx)+1;
            else % should not happen
                keyboard()
            end
        end
        if ~all(counter>1)
            for ccid=reshape(find(counter==1),1,[])
                tkey=sprintf('%d-',per_trial_motif_freq(tt,1),gnodes(conncomp==ccid));
                if ~ismember(tkey,processed)
                    disconnected=[disconnected;{per_trial_motif_freq(tt,1)},{gnodes(conncomp==ccid)}];
                    processed=[processed;tkey];
                end
            end
        end
    end
end
save(fullfile('bzdata','disconnected_motifs.mat'),"disconnected","blame")





