% tag spikes in neuron with loop or chain activity
% streamline related/duplicate code sniplet.

function run_length=chain_loop_stats(chain_replay,ring_replay,disconnected,opt)
arguments
    chain_replay = []
    ring_replay = []
    disconnected = []
    opt.skipfile (1,1) logical = true
end

if isempty(ring_replay)
    load(fullfile('binary','motif_replay.mat'),'ring_replay');
end

if isempty(chain_replay)
    load(fullfile('binary','motif_replay.mat'),'chain_replay');
end

if isempty(disconnected)
    load(fullfile('binary','nested_loops_stats.mat'),"disconnected")
end


per_spk_tag=true;

disconnKey=cell(0);
for ii=1:size(disconnected,1)
    onekey=sprintf('%d-',disconnected{ii,1},disconnected{ii,2});
    disconnKey=[disconnKey;onekey];
end


if true%~exist('inited','var') || ~inited  % denovo data generation
    inited=true;
    %% single spike chain
    ssc_sess=unique(chain_replay.session);

    %% single spike loop
    %     load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats'); % bz.rings.rings_time_constant
    ring_replay=ring_replay(ring_replay.wave~="none",:);
    ssl_sess=unique(ring_replay.session);
    usess=union(ssc_sess,ssl_sess);

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

    for sessid=double(reshape(usess,1,[]))
        len_ms=wave.replay.sessid2length(sessid)./3;
        covered=false(ceil(len_ms),1); % 2hrs of milli-sec bins
        [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
        FT_SPIKE.lc_tag=cell(size(FT_SPIKE.timestamp));

        for cidx=1:numel(FT_SPIKE.timestamp)
            FT_SPIKE.lc_tag{cidx}=zeros(size(FT_SPIKE.timestamp{cidx}),'uint8');
        end
        wtsel=trials(:,9)>0 & trials(:,10)>0 & ismember(trials(:,5),[4 8]) &ismember(trials(:,8),[3 6]);

        per_trial_motif_freq=[per_trial_motif_freq;...
            repmat(sessid,nnz(wtsel),1),find(wtsel),trials(wtsel,5),trials(wtsel,8),...
            zeros(nnz(wtsel),4)];
        per_trial_motif_freq_perwave=[per_trial_motif_freq_perwave;...
            zeros(nnz(wtsel),12)];
        per_trial_motif_cid=[per_trial_motif_cid;cell(nnz(wtsel),2)];

        %% single spike chain

        for sidx=reshape(find(chain_replay.session==sessid),1,[])
            cids=chain_replay.meta{sidx,2};

            currkey=sprintf('%d-',sessid,sort(cids));
            if ismember(currkey,disconnKey)
                warning('skipped disconnected motif') % TODO: varify logic
                continue
            end

            % ts_id: ts, cid, pos, trial_time, trial

            % per-trial frequency.
            tsel=[];
            if chain_replay.delay(sidx)==3
                switch chain_replay.wave(sidx)
                    case "s1d3"
                        tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==3;
                        perwaveidx=[2,5];
                        samp=4;
                    case "s2d3"
                        tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==3;
                        perwaveidx=[2,5];
                        samp=8;
                    case "olf_s1"
                        tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==3;
                        perwaveidx=[1,4];
                    case "olf_s2"
                        tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==3;
                        perwaveidx=[1,4];
                    case "dur_d3"
                        tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,3),[4 8]) & per_trial_motif_freq(:,4)==3;
                        perwaveidx=[3,6];
                end
            else
                switch chain_replay.wave(sidx)
                    case "s1d6"
                        tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==6;
                        perwaveidx=[2,5];
                        samp=4;
                    case "s2d6"
                        tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==6;
                        perwaveidx=[2,5];
                        samp=8;
                    case "olf_s1"
                        tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==6;
                        perwaveidx=[1,4];
                    case "olf_s2"
                        tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==6;
                        perwaveidx=[1,4];
                    case "dur_d6"
                        tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,3),[4 8]) & per_trial_motif_freq(:,4)==6;
                        perwaveidx=[3,6];
                end
            end
            if any(tsel)
                % chain motif ++
                per_trial_motif_freq(tsel,5)=per_trial_motif_freq(tsel,5)+1;
                per_trial_motif_freq_perwave(tsel,perwaveidx(1))=per_trial_motif_freq_perwave(tsel,perwaveidx(1))+1;

                for ttt=reshape(find(tsel),1,[])
                    per_trial_motif_cid{ttt,1}=[per_trial_motif_cid{ttt,1},cids];
                end
                
                trl_align=chain_replay.trl_align{sidx};
                pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                for onetid=reshape(trl_align(pref_delay,1),1,[])
                    activitysel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,2)==onetid;
                    per_trial_motif_freq(activitysel,6)=per_trial_motif_freq(activitysel,6)+1;
                    per_trial_motif_freq_perwave(activitysel,perwaveidx(2))=per_trial_motif_freq_perwave(activitysel,perwaveidx(2))+1;
                end
            end
            % ==================
            %check cid in largest network

            if per_spk_tag
                for cidx=1:size(chain_replay.ts{sidx},2)
                    cid=cids(cidx);
                    %                     ucid=[ucid,cid];
                    cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                    totag=chain_replay.ts{sidx}(:,cidx);
                    [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                    FT_SPIKE.lc_tag{cidsel}(totagidx(ism))=bitor(FT_SPIKE.lc_tag{cidsel}(totagidx(ism)),1,'uint8');
                end
            end
            % run length tag
            for cidx=1:size(chain_replay.ts{sidx},1)
                onset=floor(chain_replay.ts{sidx}(cidx,1)./3); % 0.1ms precision
                offset=ceil(chain_replay.ts{sidx}(cidx,end)./3); % 0.1ms precision
                covered(onset:offset)=true;
            end
        end


        %% single spike loop
        
        for cc=reshape(find(ring_replay.session==sessid),1,[])
            currkey=sprintf('%d-',sessid,sort(ring_replay.meta{cc,2}));
            if ismember(currkey,disconnKey)
                warning('skipped disconnected motif')
                continue
            end

            % .ts_id(:,6) => loop tag
            % per-trial frequency.
            [pref3,pref6]=bz.rings.preferred_trials_rings(ring_replay.wave(cc),trials);
            tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,2),[pref3;pref6]);
            per_trial_motif_freq(tsel,7)=per_trial_motif_freq(tsel,7)+1;
            % per wave
            % switch oneloop.rstats{5}
            %     case 'olf'
                    % perwaveidx=[7,10];
                % case 'both'
                    perwaveidx=[8,11];
                % case 'dur'
            %         perwaveidx=[9,12];
            %     otherwise
            %         keyboard();
            % end

            per_trial_motif_freq_perwave(tsel,perwaveidx(1))=per_trial_motif_freq_perwave(tsel,perwaveidx(1))+1;

            for ttt=reshape(find(tsel),1,[])
                per_trial_motif_cid{ttt,2}=[per_trial_motif_cid{ttt,2},ring_replay.meta(cc,2)];
            end
            
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);

            % TODO: WIP0828
            u_act_per_trl=unique(oneloop.ts_id(oneloop.ts_id.Loop_tag>0,{'Trial','Loop_tag'}),'rows');
            [gc,gr]=groupcounts(u_act_per_trl.Trial);
            for aa=1:numel(gr)
                asel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,2)==gr(aa);
                per_trial_motif_freq(asel,8)=per_trial_motif_freq(asel,8)+gc(aa);
                per_trial_motif_freq_perwave(asel,perwaveidx(2))=per_trial_motif_freq_perwave(asel,perwaveidx(2))+gc(aa);
            end

            %check cid in largest network

            if per_spk_tag
                for cidx=1:size(ring_replay.meta(cc,2),2)
                    % cid=ssloop_trl.meta(cc,2)(cidx);
                    %             ucid=[ucid,cid];
                    cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                    
                    totag=oneloop.ts_id.TS(oneloop.ts_id.Loop_tag>0 & oneloop.ts_id.CID==cid);
                    [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                    FT_SPIKE.lc_tag{cidsel}(totagidx)=bitor(FT_SPIKE.lc_tag{cidsel}(totagidx),4,'uint8');
                end
            end
            % run length tag
            for tagi=reshape(setdiff(unique(oneloop.ts_id.Loop_tag),0),1,[])
                tseq=oneloop.ts_id.TS(oneloop.ts_id.Loop_tag==tagi);
                onset=floor(tseq(1)./3); % 0.1ms precision
                offset=ceil(tseq(end)./3); % 0.1ms precision
                %                 if offset-onset<2,keyboard();end
                covered(onset:offset)=true;
            end
        end


        if ~opt.skipfile
            blame=vcs.blame();
            save(fullfile("bzdata","SingleSpikeChainedLoop"+num2str(sessid)+".mat"),'FT_SPIKE','covered','blame')
        end
        per_sess_coverage.("SSS"+num2str(sessid))=covered;

    end
    denovoflag=true;
else % load from file
    per_sess_coverage=struct();
    for bi=[150,300,600]
        for sessid=[14,18,22,33,34,68,100,102,114]
            load(fullfile("bzdata","SingleSpikeChainedLoop"+num2str(sessid)+".mat"),'covered')
            per_sess_coverage.("SSS"+num2str(sessid))=covered;
        end
    end
    denovoflag=false;
end
if denovoflag
    disp("New set of data file generated")
    %     keyboard();
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


%%
relax_cosec=false;
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


function relaxed_consecutive()
%%
covered=covcell(contains(covfn,"SSS"));
figure()
hold on;
cmap=colormap('lines');
cidx=1;
for consecthresh=[300 600 1500]
    run_length=[];
    for jj=1:numel(covered)
        %pass 1
        edges = find(diff([0;covered{jj};0]==1));
        onset = edges(1:2:end-1);  % Start indices
        offset = edges(2:2:end);

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
        run_length =[run_length; (offset-onset)./10];  % Consecutive ones counts
    end
    chained_loops_pdf=histcounts(run_length,[0:19,20:20:180,200:100:2000],'Normalization','pdf');
    plot([0.5:19.5,30:20:190,250:100:1950],chained_loops_pdf,'-','Color',cmap(cidx,:));

    qtrs=prctile(run_length,[10,50,90]);
    xline(qtrs,'--',string(qtrs),'Color',cmap(cidx,:))
    cidx=cidx+1;
end

xlim([2,2000])

ylim([8e-6,0.1])
set(gca(),'XScale','log','YScale','log')
xlabel('Time (ms)')
ylabel('Probability density')


end


function plot_motif_freq()
load(fullfile("bzdata","per_trial_motif_spk_freq.mat"),'per_trial_motif_freq','per_trial_motif_freq_perwave','per_trial_motif_cid')

olfsel=sum(per_trial_motif_freq_perwave(:,[1 7]),2)>0;
bothsel=sum(per_trial_motif_freq_perwave(:,[2 8]),2)>0;
olffreq=sum(per_trial_motif_freq_perwave(olfsel,[4 10]),2)./per_trial_motif_freq(olfsel,4);
bothfreq=sum(per_trial_motif_freq_perwave(bothsel,[5 11]),2)./per_trial_motif_freq(bothsel,4);
olfiqrs=prctile(olffreq,[25,50,75]);
bothiqrs=prctile(bothfreq,[25,50,75]);
fh=figure();
fh.Position(1:2)=[100,100];
hold on
boxplot([olffreq;bothfreq],[ones(size(olffreq));2*ones(size(bothfreq))],'Whisker',realmax,'Colors','k')

olfiqr=prctile(olffreq,[25,50,75,100]);
bothiqr=prctile(bothfreq,[25,50,75,100]);

xlim([0.5,2.5])
ylim([0,30])
title({sprintf('%0.2f,',olfiqr),sprintf('%0.2f,',bothiqr)});
ylabel('Motif frequency  (Hz)')
set(gca,'XTick',1:2,'XTickLabel',{'Olf','Both'})
end



%% ========================SHOWCASE===============
% moved to chain_loops_SC_spk.m

function db_test()
dbfile=fullfile("bzdata","rings_wave_burst_iter_600.db");
conn=sqlite(dbfile,"readonly");
keys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
keys(1)
conn.fetch("SELECT COUNT(*) FROM "+keys(1));
%         rids=table2array(conn.fetch("SELECT DISTINCT Var1 from d6s1d6s22r3n274_ts"));
end

