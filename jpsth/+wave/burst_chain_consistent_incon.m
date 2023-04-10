% comparison of consistent / inconsistent burst chains
% motif frequency & time constant?


burstinterval=300;
per_spk_tag=true;
bburst=true;

skipfile=false;

%% burst spike chain
sess=[];
for fprefix=["","rev_"]
    if fprefix==""
        fkey="FWD";
    else
        fkey="REV";
    end
    cfstr.(fkey)=load(fullfile("bzdata","chain_sust_"+fprefix+"tag_"+burstinterval+".mat"),'out');
    keys=[struct2cell(structfun(@(x) fieldnames(x), cfstr.(fkey).out.d6, 'UniformOutput', false));...
        struct2cell(structfun(@(x) fieldnames(x), cfstr.(fkey).out.d3, 'UniformOutput', false))];
    keys=vertcat(keys{:});
    sess=[sess;{str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once'))}];
end
usess=intersect(unique(sess{1}),unique(sess{2}));

% single spk chn:1, burst spk chn:2, single spk loop:4, burst spk loop:8

%% per-session entrance
per_sess_coverage=struct();

% 1:session#, 2:trial#, 3:sample, 4:duration, 5:chain motif count,
% 6: chain activity c6ount, 7: loop motif count, 8:loop activity count
per_trial_motif_freq=[];
% 1:olf_chain motif count, 2:both_chain motif count, 3:dur_chain motif count
% 4:olf_chain activity count, 5:both_chain activity count, 6:dur_chain activity count
% 7-12: similar, but for loops
per_trial_motif_freq_perwave=[];
per_trial_motif_cid=cell(0);

for sessid=reshape(usess,1,[])

    covered=false(1000,1); % 2hrs of milli-sec bins

    [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
    FT_SPIKE.lc_tag=cell(size(FT_SPIKE.timestamp));
    for cidx=1:numel(FT_SPIKE.timestamp)
        FT_SPIKE.lc_tag{cidx}=zeros(size(FT_SPIKE.timestamp{cidx}));
    end
    wtsel=trials(:,9)>0 & trials(:,10)>0 & ismember(trials(:,5),[4 8]) &ismember(trials(:,8),[3 6]);

    per_trial_motif_freq=[per_trial_motif_freq;...
        repmat(sessid,nnz(wtsel),1),find(wtsel),trials(wtsel,5),trials(wtsel,8),...
        zeros(nnz(wtsel),4)];
    per_trial_motif_freq_perwave=[per_trial_motif_freq_perwave;...
        zeros(nnz(wtsel),12)];
    per_trial_motif_cid=[per_trial_motif_cid;cell(nnz(wtsel),2)];


    %% multi spike chain
    fkey="FWD";
    for dur=["d6","d3"]
        for wid=reshape(fieldnames(cfstr.(fkey).out.(dur)),1,[])
            for cc=reshape(fieldnames(cfstr.(fkey).out.(dur).(wid{1})),1,[])
                if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                    continue
                end
                tsel=[];
                if dd==3
                    switch wid{1}
                        case 's1d3'
                            tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==3;
                            perwaveidx=[2,5];
                        case 's2d3'
                            tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==3;
                            perwaveidx=[2,5];
                        case 'olf_s1'
                            tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==3;
                            perwaveidx=[1,4];
                        case 'olf_s2'
                            tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==3;
                            perwaveidx=[1,4];
                        case 'dur_d3'
                            tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,3),[4 8]) & per_trial_motif_freq(:,4)==3;
                            perwaveidx=[3,6];
                    end
                else
                    switch wid{1}
                        case 's1d6'
                            tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==6;
                            perwaveidx=[2,5];
                        case 's2d6'
                            tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==6;
                            perwaveidx=[2,5];
                        case 'olf_s1'
                            tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==6;
                            perwaveidx=[1,4];
                        case 'olf_s2'
                            tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==6;
                            perwaveidx=[1,4];
                        case 'dur_d6'
                            tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,3),[4 8]) & per_trial_motif_freq(:,4)==6;
                            perwaveidx=[3,6];
                    end
                end
                if ~isempty(tsel)
                    % chain motif ++
                    per_trial_motif_freq(tsel,5)=per_trial_motif_freq(tsel,5)+1;
                    per_trial_motif_freq_perwave(tsel,perwaveidx(1))=per_trial_motif_freq_perwave(tsel,perwaveidx(1))+1;

                    for ttt=reshape(find(tsel),1,[])
                        per_trial_motif_cid{ttt,1}=[per_trial_motif_cid{ttt,1},onechain.meta(1)];
                    end

                    [~,tid]=ismember(onechain.ts(:,1),onechain.ts_id(:,1));
                    for onetid=reshape(tid,1,[])
                        activitysel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,2)==onechain.ts_id(onetid,5);
                        per_trial_motif_freq(activitysel,6)=per_trial_motif_freq(activitysel,6)+1;
                        per_trial_motif_freq_perwave(activitysel,perwaveidx(2))=per_trial_motif_freq_perwave(activitysel,perwaveidx(2))+1;
                    end
                end

                if per_spk_tag
                    for cidx=1:size(onechain.meta{1},2)
                        cid=onechain.meta{1}(cidx);
                        %                     ucid=[ucid,cid];
                        cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                        for bscid=1:numel(onechain.ts)
                            totag=onechain.ts{bscid}(onechain.ts{bscid}(:,1)==cidx,3);
                            [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                            FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+2;
                        end
                    end
                end
                % run length tag
                for bscid=1:numel(onechain.ts)
                    onset=floor(onechain.ts{bscid}(1,3)./3); % 0.1ms precision
                    offset=ceil(onechain.ts{bscid}(end,3)./3); % 0.1ms precision
                    covered(onset:offset)=true;
                end
            end
        end
    end





    if ~skipfile
        blame=vcs.blame();
        if bburst
            save(fullfile("bzdata","ChainedLoop"+burstinterval+"S"+num2str(sessid)+".mat"),'FT_SPIKE','covered','blame')
        else
            save(fullfile("bzdata","SingleSpikeChainedLoop"+num2str(sessid)+".mat"),'FT_SPIKE','covered','blame')
        end
    end
    
    per_sess_coverage.("B"+burstinterval+"S"+num2str(sessid))=covered;

end
denovoflag=true;
% else % load from file
%     per_sess_coverage=struct();
%     for bi=[150,300,600]
%         for sessid=[14,18,22,33,34,68,100,102,114]
%             if bburst
%                 load(fullfile("bzdata","ChainedLoop"+bi+"S"+num2str(sessid)+".mat"),'covered','FT_SPIKE')
%             else
%                 load(fullfile("bzdata","SingleSpikeChainedLoop"+num2str(sessid)+".mat"),'covered')
%             end
%             per_sess_coverage.("B"+bi+"S"+num2str(sessid))=covered;
%         end
%     end
%     denovoflag=false;

if denovoflag
    disp("New set of data file generated")
    keyboard();
    if false
        blame=vcs.blame();
        save(fullfile("bzdata","per_trial_motif_spk_freq.mat"),'per_trial_motif_freq','per_trial_motif_freq_perwave','per_trial_motif_cid','blame')
        % chains
        unique(cellfun(@(x) numel(x),per_trial_motif_cid(:,1)))
        % loops
        trlcnt=cellfun(@(x) numel(x),per_trial_motif_cid(:,2))
        unique(trlcnt(per_trial_motif_freq(:,7)>1))
        % chained-loops
        trlcnt=arrayfun(@(x) numel(unique([per_trial_motif_cid{x,1},per_trial_motif_cid{x,2}])),1:size(per_trial_motif_cid,1))
        unique(trlcnt(per_trial_motif_freq(:,5)>=1 & per_trial_motif_freq(:,7)>=1))

    end
end
covcell=struct2cell(per_sess_coverage);
covfn=fieldnames(per_sess_coverage);
covered=struct();
run_length=struct();

if bburst
    %%
    for bi=600%[150,300,600]
        covered.("B"+bi)=covcell(contains(covfn,"B"+bi));
        run_length.("B"+bi)=[];
        for jj=1:numel(covered.("B"+bi))
            edges = find(diff([0;covered.("B"+bi){jj};0]==1));
            onset = edges(1:2:end-1);  % Start indices
            %     disp([jj,min(edges(2:2:end)-onset)]);
            run_length.("B"+bi) =[run_length.("B"+bi); edges(2:2:end)-onset];  % Consecutive ones counts
            %     per_sess_coverage.("S"+num2str(sessid))=run_length;
        end
    end
    figure()
    hold on;
    ph=[];
    cmap=colormap("lines");
    cidx=1;
    for bi=[150,300,600]
        chained_loops_pdf=histcounts(run_length.("B"+bi),[0:2:18,20:20:100,200,300:300:1200],'Normalization','pdf');
        ph(end+1)=plot([1:2:19,30:20:90,150,250,450:300:1050],chained_loops_pdf,'-','Color',cmap(cidx,:));
        qtrs=prctile(run_length.("B"+bi),[25,50,75]);
        %     qtrs19=prctile(run_length,[10,20,80,90]);
        xline(qtrs,'--',string(qtrs),'Color',cmap(cidx,:)) % 17 24 35
        cidx=cidx+1;
    end
    xlim([5,1200])
    ylim([8e-7,0.1])
    set(gca(),'XScale','log','YScale','log')
    xlabel('Time (ms)')
    ylabel('Probability density')
else
    %%
    relax_cosec=true;
    consecthresh=1500;
    covered=covcell(contains(covfn,"SSS"));
    run_length=[];
    for jj=1:numel(covered)
        %pass 1
        edges = find(diff([0;covered{jj};0]==1));
        onset = edges(1:2:end-1);  % Start indices
        offset = edges(2:2:end);

        if relax_cosec
            % patch through
            latmat=onset-offset.'; % row->onset, col->offset
            for onidx=1:numel(onset)
                consec=find(latmat(onidx,:)>0 & latmat(onidx,:)<=consecthresh);
                if ~isempty(consec)
                    offidx=min(consec);
                    covered{jj}((offset(offidx)-1):onset(onidx))=1;
                end
            end
            % pass 2
            edges = find(diff([0;covered{jj};0]==1));
            onset = edges(1:2:end-1);  % Start indices
            offset = edges(2:2:end);
        end
        run_length =[run_length; (offset-onset)./10];  % Consecutive ones counts
    end

    figure()
    hold on;
    if relax_cosec
        chained_loops_pdf=histcounts(run_length,[0:19,20:20:180,200:100:2000],'Normalization','pdf');
        plot([0.5:19.5,30:20:190,250:100:1950],chained_loops_pdf,'-k');
        xlim([2,2000])
        titie("Consecutive window "+(consecthresh/30)+" msec")
    else
        chained_loops_pdf=histcounts(run_length,[0:19,20:20:300],'Normalization','pdf');
        plot([0.5:19.5,30:20:290],chained_loops_pdf,'-k');
        xlim([2,500])
    end
    qtrs=prctile(run_length,[10,50,90]);
    xline(qtrs,'--k',string(qtrs)) % 17 24 35

    ylim([8e-6,0.1])
    set(gca(),'XScale','log','YScale','log')
    xlabel('Time (ms)')
    ylabel('Probability density')
end
