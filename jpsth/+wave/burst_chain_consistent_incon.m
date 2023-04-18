% comparison of consistent / inconsistent burst chains
% motif frequency
% spike proportion, i.e. revert spike tag pipeline
% time constant
skipplot=true;
appear=struct();

for burstinterval=[150 300 600]
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

    %% per-session entrance
    per_sess_coverage=struct();

    % 1:session#, 2:trial#, 3:sample, 4:duration, 5:chain motif count,
    % 6: chain activity c6ount, 7: loop motif count, 8:loop activity count
    per_trial_motif_freq=cell2struct({[],[]},{'FWD','REV'},2);
    % 1:olf_chain motif count, 2:both_chain motif count, 3:dur_chain motif count
    % 4:olf_chain activity count, 5:both_chain activity count, 6:dur_chain activity count
    % 7-12: similar, but for loops
    per_trial_motif_freq_perwave=cell2struct({[],[]},{'FWD','REV'},2);
    per_trial_motif_cid=cell2struct({cell(0),cell(0)},{'FWD','REV'},2);

    lc_tags=cell2struct({[],[]},{'FWD','REV'},2);

    for sessid=reshape(usess,1,[])

        [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
        sess_lc_tag=cell(size(FT_SPIKE.timestamp));
        for cidx=1:numel(FT_SPIKE.timestamp)
            sess_lc_tag{cidx}=zeros(size(FT_SPIKE.timestamp{cidx}),'uint8');
        end

        wtsel=trials(:,9)>0 & trials(:,10)>0 & ismember(trials(:,5),[4 8]) &ismember(trials(:,8),[3 6]);

        %% multi spike chain

        for fkey=["FWD","REV"]
            per_trial_motif_freq.(fkey)=[per_trial_motif_freq.(fkey);...
                repmat(sessid,nnz(wtsel),1),find(wtsel),trials(wtsel,5),trials(wtsel,8),...
                zeros(nnz(wtsel),4)];
            per_trial_motif_freq_perwave.(fkey)=[per_trial_motif_freq_perwave.(fkey);...
                zeros(nnz(wtsel),12)];
            per_trial_motif_cid.(fkey)=[per_trial_motif_cid.(fkey);cell(nnz(wtsel),2)];

            covered=false(1000,1); % 2hrs of milli-sec bins
            for dur=["d6","d3"]
                for wid=reshape(fieldnames(cfstr.(fkey).out.(dur)),1,[])
                    for cc=reshape(fieldnames(cfstr.(fkey).out.(dur).(wid{1})),1,[])
                        if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                            continue
                        end
                        tsel=[];
                        if dur=="d3"
                            switch wid{1}
                                case 's1d3'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,3)==4 & per_trial_motif_freq.(fkey)(:,4)==3;
                                    perwaveidx=[2,5];
                                case 's2d3'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,3)==8 & per_trial_motif_freq.(fkey)(:,4)==3;
                                    perwaveidx=[2,5];
                                case 'olf_s1'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,3)==4 & per_trial_motif_freq.(fkey)(:,4)==3;
                                    perwaveidx=[1,4];
                                case 'olf_s2'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,3)==8 & per_trial_motif_freq.(fkey)(:,4)==3;
                                    perwaveidx=[1,4];
                                case 'dur_d3'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & ismember(per_trial_motif_freq.(fkey)(:,3),[4 8]) & per_trial_motif_freq.(fkey)(:,4)==3;
                                    perwaveidx=[3,6];
                            end
                        else
                            switch wid{1}
                                case 's1d6'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,3)==4 & per_trial_motif_freq.(fkey)(:,4)==6;
                                    perwaveidx=[2,5];
                                case 's2d6'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,3)==8 & per_trial_motif_freq.(fkey)(:,4)==6;
                                    perwaveidx=[2,5];
                                case 'olf_s1'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,3)==4 & per_trial_motif_freq.(fkey)(:,4)==6;
                                    perwaveidx=[1,4];
                                case 'olf_s2'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,3)==8 & per_trial_motif_freq.(fkey)(:,4)==6;
                                    perwaveidx=[1,4];
                                case 'dur_d6'
                                    tsel=per_trial_motif_freq.(fkey)(:,1)==sessid & ismember(per_trial_motif_freq.(fkey)(:,3),[4 8]) & per_trial_motif_freq.(fkey)(:,4)==6;
                                    perwaveidx=[3,6];
                            end
                        end
                        if ~isempty(tsel)
                            onechain=cfstr.(fkey).out.(dur).(wid{1}).(cc{1});

                            % chain motif ++
                            per_trial_motif_freq.(fkey)(tsel,5)=per_trial_motif_freq.(fkey)(tsel,5)+1;
                            per_trial_motif_freq_perwave.(fkey)(tsel,perwaveidx(1))=per_trial_motif_freq_perwave.(fkey)(tsel,perwaveidx(1))+1;

                            for ttt=reshape(find(tsel),1,[])
                                per_trial_motif_cid.(fkey){ttt,1}=[per_trial_motif_cid.(fkey){ttt,1},onechain.meta(1)];
                            end

                            [~,tid]=ismember(cellfun(@(x) x(1,3), onechain.ts),onechain.ts_id(:,1));
                            for onetid=reshape(tid,1,[])
                                activitysel=per_trial_motif_freq.(fkey)(:,1)==sessid & per_trial_motif_freq.(fkey)(:,2)==onechain.ts_id(onetid,5);
                                per_trial_motif_freq.(fkey)(activitysel,6)=per_trial_motif_freq.(fkey)(activitysel,6)+1;
                                per_trial_motif_freq_perwave.(fkey)(activitysel,perwaveidx(2))=per_trial_motif_freq_perwave.(fkey)(activitysel,perwaveidx(2))+1;
                            end
                        end
                        % per_spk_tag

                        for cidx=1:size(onechain.meta{1},2)
                            cid=onechain.meta{1}(cidx);
                            %                     ucid=[ucid,cid];
                            cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                            for bscid=1:numel(onechain.ts)
                                burstsel=~([diff(onechain.ts{bscid}(:,1));1]);
                                burst_tag=onechain.ts{bscid}(onechain.ts{bscid}(:,1)==cidx & burstsel,3);
                                [~,burst_tagidx]=ismember(burst_tag,FT_SPIKE.timestamp{cidsel});
                                sess_lc_tag{cidsel}(burst_tagidx)=bitor(sess_lc_tag{cidsel}(burst_tagidx),2,"uint8");
                                fc_tag=onechain.ts{bscid}(onechain.ts{bscid}(:,1)==cidx & ~burstsel,3);
                                [~,fc_tagidx]=ismember(fc_tag,FT_SPIKE.timestamp{cidsel});
                                sess_lc_tag{cidsel}(fc_tagidx)=bitor(sess_lc_tag{cidsel}(fc_tagidx),1,"uint8");
                            end
                        end

                        % run_length
                        for bscid=1:numel(onechain.ts)
                            onset=floor(onechain.ts{bscid}(1,3)./3); % 0.1ms precision
                            offset=ceil(onechain.ts{bscid}(end,3)./3); % 0.1ms precision
                            covered(onset:offset)=true;
                        end
                    end
                end
            end
            % could save per-session data here
            % considered unnecessary
            per_sess_coverage.(fkey).("B"+burstinterval+"S"+num2str(sessid))=covered;
            lc_tags.(fkey).("S"+sessid)=sess_lc_tag;
        end
    end


    %% appearance
    ukeys=cell2struct({[],[]},{'FWD','REV'},2);
    for fkey=["FWD","REV"]
        for ii=1:size(per_trial_motif_freq.(fkey),1)
            if ~isempty(per_trial_motif_cid.(fkey){ii,1})
                ukeys.(fkey)=[ukeys.(fkey),cellfun(@(x) string(sprintf('%d-',per_trial_motif_freq.(fkey)(ii,1),x)),per_trial_motif_cid.(fkey){ii,1})];
            end
        end
    end
    ukeys.UNIQFWD=unique(ukeys.FWD);
    ukeys.UNIQREV=unique(ukeys.REV);
    appear.("B"+burstinterval).FWD=numel(ukeys.UNIQFWD);
    appear.("B"+burstinterval).REV=numel(ukeys.UNIQREV);

    %% per_chainid_repeats
    %     fsel=per_trial_motif_freq.FWD(:,5)>0;
    %     mean(per_trial_motif_freq.FWD(fsel,6)./per_trial_motif_freq.FWD(fsel,5))
    %     rsel=per_trial_motif_freq.REV(:,5)>0;
    %     mean(per_trial_motif_freq.REV(rsel,6)./per_trial_motif_freq.REV(rsel,5))
    if ~skipplot
        figure();
        tiledlayout('flow')
        nexttile
        bar([numel(ukeys.UNIQFWD),numel(ukeys.UNIQREV)])
        ylabel('Number of chains with burst sequence')
        set(gca,'XTick',1:2,'XTickLabel',{'Consistent','Inconsistent'})
        title('Unique chain with spike activity')

        %% time constant
        %     figure()
        nexttile
        hold on;
        ph=[];
        cmap=[1,0,0;0,0,0];
        cidx=1;

        run_length=struct();
        covered=struct();

        for fkey=["FWD","REV"]
            covcell=struct2cell(per_sess_coverage.(fkey));
            covfn=fieldnames(per_sess_coverage.(fkey));

            covered.("B"+burstinterval).(fkey)=covcell(contains(covfn,"B"+burstinterval));
            run_length.("B"+burstinterval).(fkey)=[];
            for jj=1:numel(covered.("B"+burstinterval).(fkey))
                edges = find(diff([0;covered.("B"+burstinterval).(fkey){jj};0]==1));
                onset = edges(1:2:end-1);  % Start indices
                %     disp([jj,min(edges(2:2:end)-onset)]);
                run_length.("B"+burstinterval).(fkey) =[run_length.("B"+burstinterval).(fkey); (edges(2:2:end)-onset)./10];  % Consecutive ones counts
                %     per_sess_coverage.("S"+num2str(sessid))=run_length;
            end

            chained_loops_pdf=histcounts(run_length.("B"+burstinterval).(fkey),[0:2:18,20:20:100,200,300:300:1200],'Normalization','pdf');
            ph(cidx)=plot([1:2:19,30:20:90,150,250,450:300:1050],chained_loops_pdf,'-','Color',cmap(cidx,:));
            qtrs=prctile(run_length.("B"+burstinterval).(fkey),[10,50,90]);
            %     qtrs19=prctile(run_length,[10,20,80,90]);
            xline(qtrs,'--',string(qtrs),'Color',cmap(cidx,:)) % 17 24 35
            cidx=cidx+1;
        end
        p=ranksum(run_length.("B"+burstinterval).FWD,run_length.("B"+burstinterval).REV);
        title(sprintf('p = %0.3f',p))
        legend(ph,{'Consistent','Inconsistent'})
        xlim([5,1200])
        ylim([8e-7,0.1])
        set(gca(),'XScale','log','YScale','log')
        xlabel('Time (ms)')
        ylabel('Probability density')


        %% motif freq

        % function plot_motif_freq()
        % load(fullfile("bzdata","per_trial_motif_spk_freq.mat"),'per_trial_motif_freq','per_trial_motif_freq_perwave','per_trial_motif_cid')

        %     fh=figure();
        %     fh.Position(1:2)=[100,100];
        %     tiledlayout(1,2)
        aax=struct();
        for fkey=["FWD","REV"]
            olfsel=per_trial_motif_freq_perwave.(fkey)(:,1)>0;
            bothsel=per_trial_motif_freq_perwave.(fkey)(:,2)>0;
            olffreq.(fkey)=per_trial_motif_freq_perwave.(fkey)(olfsel,4)./per_trial_motif_freq.(fkey)(olfsel,4); % col4=>duration
            bothfreq.(fkey)=per_trial_motif_freq_perwave.(fkey)(bothsel,5)./per_trial_motif_freq.(fkey)(bothsel,4);
            %         continue
            aax.(fkey)=nexttile();
            hold on
            %         olfiqrs=prctile(olffreq.(fkey),[25,50,75]);
            %         bothiqrs=prctile(bothfreq.(fkey),[25,50,75]);

            boxplot([olffreq.(fkey);bothfreq.(fkey)],[ones(size(olffreq.(fkey)));2*ones(size(bothfreq.(fkey)))],'Colors','k','OutlierSize',4,'Symbol','c+')

            olfiqr=prctile(olffreq.(fkey),[25,50,75,100]);
            bothiqr=prctile(bothfreq.(fkey),[25,50,75,100]);

            spanolf=olfiqr(3)+1.5*(diff(olfiqr(1:2:3)));
            spanboth=bothiqr(3)+1.5*(diff(bothiqr(1:2:3)));

            yspan.(fkey)=max([spanboth;spanolf])*1.25;

            xlim([0.5,2.5])
            ylabel('Motif frequency  (Hz)')
            set(gca,'XTick',1:2,'XTickLabel',{'Olf','Both'})
            if fkey=="FWD"
                title("Total consistent motif freq");
            else
                title("Total inconsistent motif freq");
            end
        end

        cellfun(@(x) set(x,'YLim',[0,max(struct2array(yspan))]),struct2cell(aax))



        %% spike proportion

        sums=cell2struct({[],[]},{'FWD','REV'},2);
        for fkey=["FWD","REV"]
            sessfn=fieldnames(lc_tags.(fkey));
            for fn=reshape(sessfn,1,[])
                sessid=str2double(replace(fn{1},'S',''));
                [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);

                sesssel=per_trial_motif_freq.(fkey)(:,1)==sessid;
                sess_uid=unique(cell2mat([per_trial_motif_cid.(fkey){sesssel,1}]));
                for onecid=sess_uid
                    % trial sel
                    trialsel=cellfun(@(x) ismember(onecid,cell2mat(x)), per_trial_motif_cid.(fkey)(sesssel,1));
                    trial3sel=intersect(find(trialsel), find(FT_SPIKE.trialinfo(:,8)==3));
                    trial6sel=intersect(find(trialsel), find(FT_SPIKE.trialinfo(:,8)==6));

                    idsel=strcmp(FT_SPIKE.label,num2str(onecid));
                    tssel=(ismember(FT_SPIKE.trial{idsel},trial3sel) & FT_SPIKE.time{idsel}>1 & FT_SPIKE.time{idsel}<=4)...
                        |(ismember(FT_SPIKE.trial{idsel},trial6sel) & FT_SPIKE.time{idsel}>1 & FT_SPIKE.time{idsel}<=7);
                    %                 if any(lc_tags.(ftype).(fn{1}){idsel}(tssel)==4)
                    %                     disp(4)
                    %                     keyboard()
                    %                 end

                    fc_tagged=nnz(bitand(lc_tags.(fkey).(fn{1}){idsel}(tssel),1));
                    burst_tagged=nnz(bitand(lc_tags.(fkey).(fn{1}){idsel}(tssel),2));
                    % 1:sess, 2:cid, 3:chain fc spk count, 4: chain burst spk count, 5:total spk count
                    sums.(fkey)=[sums.(fkey);sessid,onecid,fc_tagged,burst_tagged,nnz(tssel)];
                end
            end
        end
        % mean
        %     [mean(sums.FWD(:,3:4)./sums.FWD(:,5)),mean(sums.REV(:,3:4)./sums.REV(:,5))]
        %     % median
        %     [median(sums.FWD(:,3:4)./sums.FWD(:,5)),median(sums.REV(:,3:4)./sums.REV(:,5))]

        boxdata=[sums.FWD(:,3)./sums.FWD(:,5),ones(size(sums.FWD(:,3)));sums.REV(:,3)./sums.REV(:,5),2*ones(size(sums.REV(:,3)))];
        nexttile()
        boxplot(boxdata(:,1),boxdata(:,2),'Colors','k','Symbol','c+','OutlierSize',4)
        %     ylim([0,0.05])

        prcdata=[prctile(boxdata(boxdata(:,2)==1),[25,75]);...
            prctile(boxdata(boxdata(:,2)==2),[25,75])];
        ylim([0,1.25*max(prcdata(:,2)+1.5*diff(prcdata,1,2))])

        set(gca,'XTick',1:2,'XTickLabel',{'Consistent','inconsistent'},'YTick',0:0.02:0.04,'YTickLabel',0:2:4)
        ylabel('Proportion of associated spikes (%)')
        title('Proportion of spikes')


        sgtitle('Burst interval')
    end
end

figure()
bh=bar([appear.B150.FWD,appear.B300.FWD,appear.B600.FWD;appear.B150.REV,appear.B300.REV,appear.B600.REV].');
bh(1).FaceColor='w';
bh(2).FaceColor='k';
set(gca,'XTick',1:3,'XTickLabel',{'200Hz','100Hz','50Hz'})
ylabel('Number of occurrence')
legend(bh,{'Consistent','Inconsistent'},'Location','northoutside','Orientation','horizontal')
title('Burst chain occurrence')



