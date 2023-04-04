% tag spikes in neuron with loop or chain activity
per_spk_tag=true;
bburst=false;
burstinterval=600;
skipfile=true;

[sig,~]=bz.load_sig_sums_conn_file('pair',false);

optmem=false;
if ispc
    [~,sysmem]=memory();
    if sysmem.PhysicalMemory.Total < 2e10
        optmem=true;
    end
end

if true%~exist('inited','var') || ~inited  % denovo data generation
    inited=true;
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

    if bburst
        %% burst spike chain
        bschain=load(fullfile("bzdata","chain_sust_tag_"+burstinterval+".mat"),'out');
        keys=[struct2cell(structfun(@(x) fieldnames(x), bschain.out.d6, 'UniformOutput', false));...
            struct2cell(structfun(@(x) fieldnames(x), bschain.out.d3, 'UniformOutput', false))];
        keys=vertcat(keys{:});
        bsc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));

        %% burst spike loop, keys only
        dbfile=fullfile("bzdata","rings_wave_burst_iter_"+burstinterval+".db");
        conn=sqlite(dbfile,"readonly");
        bslkeys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
        close(conn);
        bsl_sess=unique(str2double(regexp(bslkeys,'(?<=s)\d{1,3}(?=r)','match','once')));

        usess=intersect(intersect(intersect(ssc_sess,bsc_sess),ssl_sess),bsl_sess);
    end
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

    for sessid=[14,18]%reshape(usess,1,[])
        % extract connected components with graph tools;
        ssel=sig.sess==sessid;
        gh=graph(cellstr(int2str(sig.suid(ssel,1))),cellstr(int2str(sig.suid(ssel,2))));
        [sgbin,binsize]=conncomp(gh);
        bidx=binsize(sgbin)==max(binsize);
        largec=subgraph(gh,bidx);
        conn_cid=str2double(table2cell(largec.Nodes));
        %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
        coveredNG=false(1000,1); % 2hrs of milli-sec bins

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
                        activitysel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,2),onechain.ts_id(tid,5));
                        per_trial_motif_freq(activitysel,6)=per_trial_motif_freq(activitysel,6)+1;
                        per_trial_motif_freq_perwave(activitysel,perwaveidx(2))=per_trial_motif_freq_perwave(activitysel,perwaveidx(2))+1;
                    end
                    % ==================
                    %check cid in largest network
                    if ~all(ismember(onechain.meta{1},conn_cid),'all')
                        disp("Outside largest component")
                        continue
                    end
                    if per_spk_tag
                        for cidx=1:size(onechain.ts,2)
                            cid=onechain.meta{1}(cidx);
                            %                     ucid=[ucid,cid];
                            cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                            totag=onechain.ts(:,cidx);
                            [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                            FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+1;
                        end
                    end
                    % run length tag
                    for cidx=1:size(onechain.ts,1)
                        onset=floor(onechain.ts(cidx,1)./30);
                        offset=ceil(onechain.ts(cidx,end)./30);
                        covered(onset:offset)=true;
                    end
                end
            end
        end
        if optmem
            clear sschain
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
            per_trial_motif_freq(tsel,7)=per_trial_motif_freq(tsel,7)+1;
            % per wave
            if strcmp(onechain.rstats{5},'olf')
                perwaveidx=[7,10];
            elseif strcmp(onechain.rstats{5},'dur')    
                perwaveidx=[9,12];

            for ttt=reshape(find(tsel),1,[])
                per_trial_motif_cid{ttt,2}=unique([per_trial_motif_cid{ttt,2},onechain.rstats{3}]);
            end

            u_act_per_trl=unique(onechain.ts_id(onechain.ts_id(:,6)>0,[5 6]),'rows');
            [gc,gr]=groupcounts(u_act_per_trl(:,1));
            for aa=1:numel(gr)
                asel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,2)==gr(aa);
                per_trial_motif_freq(asel,8)=per_trial_motif_freq(asel,8)+gc(aa);
            end
            
            %check cid in largest network
            if ~all(ismember(onechain.rstats{3},conn_cid),'all')
                disp("Outside largest component")
                continue
            end
            if per_spk_tag
                for cidx=1:size(onechain.rstats{3},2)
                    cid=onechain.rstats{3}(cidx);
                    %             ucid=[ucid,cid];
                    cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                    totag=onechain.ts_id(onechain.ts_id(:,6)>0 & onechain.ts_id(:,2)==cid,1);
                    [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                    FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+4;
                end
            end
            % run length tag
            for tagi=reshape(setdiff(unique(onechain.ts_id(:,6)),0),1,[])
                tseq=onechain.ts_id(onechain.ts_id(:,6)==tagi,1);
                onset=floor(tseq(1)./30);
                offset=ceil(tseq(end)./30);
                %                 if offset-onset<2,keyboard();end
                covered(onset:offset)=true;
            end
        end
        if optmem
            clear pstats
        end
        if bburst
            %% multi spike chain
            for dur=["d6","d3"]
                for wid=reshape(fieldnames(bschain.out.(dur)),1,[])
                    for cc=reshape(fieldnames(bschain.out.(dur).(wid{1})),1,[])
                        if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                            continue
                        end
                        %check cid in largest network
                        onechain=bschain.out.(dur).(wid{1}).(cc{1});
                        if ~all(ismember(onechain.meta{1},conn_cid),'all')
                            disp("Outside largest component")
                            continue
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
                            onset=floor(onechain.ts{bscid}(1,3)./30);
                            offset=ceil(onechain.ts{bscid}(end,3)./30);
                            covered(onset:offset)=true;
                        end
                    end
                end
            end

            if optmem
                clear bschain
            end

            %% burst spike loop
            dbfile=fullfile("bzdata","rings_wave_burst_iter_"+burstinterval+".db");
            conn=sqlite(dbfile,"readonly");
            for cc=reshape(bslkeys,1,[])
                if ~contains(cc,['s',num2str(sessid),'r']) || ~endsWith(cc,'_ts')
                    continue
                end
                chainmeta=table2array(conn.sqlread(replace(cc,'_ts','_meta')));
                %check cid in largest network
                if ~all(ismember(chainmeta,conn_cid),'all')
                    disp("Outside largest component")
                    continue
                end
                %   Not enough memory for larger complete tables
                maxrid=table2array(conn.fetch("SELECT MAX(Var1) from "+cc));
                for rid=1:1000:maxrid
                    onechain=table2array(conn.fetch("SELECT * FROM "+cc+" WHERE Var1>="+num2str(rid)+ " AND Var1<"+num2str(rid+1000)));
                    if per_spk_tag
%                         disp(rid);
                        for cidx=1:numel(chainmeta)
                            cid=chainmeta(cidx);
                            cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                            totag=onechain(onechain(:,2)==cidx,4);
                            [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                            FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+8;
                        end
                    end
                    for rrid=rid:rid+999
                        if rrid>maxrid
                            break
                        end
                        tseq=onechain(onechain(:,1)==rrid,4);
                        onset=floor(tseq(1)./30);
                        offset=ceil(tseq(end)./30);
                        covered(onset:offset)=true;
                    end
                end
            end
            conn.close();
        end
        if false % WIP % remove unconnected motifs % no effect on time constant
            underthrestrl=per_trial_motif_freq(per_trial_motif_freq(:,1)==sessid...
                & ~(per_trial_motif_freq(:,5)>=1 & per_trial_motif_freq(:,7)>=1)...
                ,2);
            for tt=reshape(underthrestrl,1,[])
                covered(floor(trials(tt,1)./30)+(0:7000))=0;
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
else % load from file
    per_sess_coverage=struct();
    for bi=[150,300,600]
        for sessid=[14,18,22,33,34,68,100,102,114]
            if bburst
                load(fullfile("bzdata","ChainedLoop"+bi+"S"+num2str(sessid)+".mat"),'covered','FT_SPIKE')
            else
                load(fullfile("bzdata","SingleSpikeChainedLoop"+num2str(sessid)+".mat"),'covered')
            end
            per_sess_coverage.("B"+bi+"S"+num2str(sessid))=covered;
        end
    end
    denovoflag=false;
end
if denovoflag
    disp("New set of data file generated")
    keyboard();
    if false
        save(fullfile("bzdata","per_trial_motif_spk_freq.mat"),'per_trial_motif_freq')
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

if bburst
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
else
    figure()
    hold on;
    chained_loops_pdf=histcounts(run_length,[0:19,20:20:300],'Normalization','pdf');
    plot([0.5:19.5,30:20:290],chained_loops_pdf,'-k');
    qtrs=prctile(run_length,[25,50,75]);
    xline(qtrs,'--k',string(qtrs)) % 17 24 35
    xlim([3,300])
    ylim([8e-6,0.1])

end
set(gca(),'XScale','log','YScale','log')
xlabel('Time (ms)')
ylabel('Probability density')


function plot_motif_freq()
load(fullfile("bzdata","per_trial_motif_spk_freq.mat"))
if false
    motif_n=per_trial_motif_freq(:,5)+per_trial_motif_freq(:,7);
    figure();histogram(motif_n)

    boxdata=[];
    sel01_100=motif_n>0 & motif_n<=100;
    boxdata=[boxdata;(per_trial_motif_freq(sel01_100,6)+per_trial_motif_freq(sel01_100,8))./per_trial_motif_freq(sel01_100,4),ones(nnz(sel01_100),1)];

    sel100_300=motif_n>100 & motif_n<=300;
    boxdata=[boxdata;(per_trial_motif_freq(sel100_300,6)+per_trial_motif_freq(sel100_300,8))./per_trial_motif_freq(sel100_300,4),ones(nnz(sel100_300),1)*2];

    sel300_500=motif_n>300 & motif_n<=500;
    boxdata=[boxdata;(per_trial_motif_freq(sel300_500,6)+per_trial_motif_freq(sel300_500,8))./per_trial_motif_freq(sel300_500,4),ones(nnz(sel300_500),1)*3];

    sel500_700=motif_n>500 & motif_n<=700;
    boxdata=[boxdata;(per_trial_motif_freq(sel500_700,6)+per_trial_motif_freq(sel500_700,8))./per_trial_motif_freq(sel500_700,4),ones(nnz(sel500_700),1)*4];

    figure('Position',[100,100,300,300])
    boxplot(boxdata(:,1),boxdata(:,2),'Colors','k','Whisker',realmax)
    set(gca(),'XTick',1:4,'XTickLabel',{'1-100','100-300','300-500','500-700'},'XTickLabelRotation',90)
    ylabel('Motif frequency (Hz)')
    xlabel('Number of motifs')
end
msel=(per_trial_motif_freq(:,5)+per_trial_motif_freq(:,7))>0;
mfreq=(per_trial_motif_freq(msel,6)+per_trial_motif_freq(msel,8))./per_trial_motif_freq(msel,4);
pct90=prctile(mfreq,[10,50,90]);
fh=figure();
tiledlayout(1,2)
nexttile
fh.Position(1:2)=[100,100];
hold on
swarmchart(ones(size(mfreq)),mfreq,'k.')
set(gca,'YScale','log')
nexttile
hold on
boxplot(mfreq,'Whisker',realmax,'Colors','k')
plot([0.9,1.1],[31,31],'k-')
xlim([0.5,1.5])
% set(gca,'YScale','log')
ylim([0,32])
title(string(pct90))
text(1,30,'max @ 221','VerticalAlignment','top','HorizontalAlignment','center')
ylabel('Motif frequency  (Hz)')
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


