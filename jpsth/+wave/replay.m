classdef replay < handle
    methods(Static)
        function [motif_trl,sum_stats,raw]=stats(motif_trl,opt)
            arguments
                motif_trl
                opt.var_len (1,1) logical = false
                % opt.nonmem_ring (1,1) logical = false
            end
            sps=30000;

            for dd=reshape(fieldnames(motif_trl),1,[])
                for ww=reshape(fieldnames(motif_trl.(dd{1})),1,[])
                    for cc=reshape(fieldnames(motif_trl.(dd{1}).(ww{1})),1,[])
                        onechain=motif_trl.(dd{1}).(ww{1}).(cc{1});
                        if ~strcmp(ww{1},'none')
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
                        else
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
                            & trl_align(:,2)>=(trl_align(:,4)+4)...  % not in delay or test
                            & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
                        freqstats.pref_succeed_ITI=[sum(pref_succeed_iti.*len),...
                            sum((onechain.trials(find(pref_delay_trls)+1,1)-onechain.trials(pref_delay_trls,2))./sps-3)]; %rwd + test

                        if ~nnonmem
                            % succeed ITI pref error
                            pref_succeed_iti_err=trl_align(:,6)==0 & trl_align(:,7)==1 ...
                                & trl_align(:,2)>=(trl_align(:,4)+4)...  % not in delay or test
                                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
                            freqstats.pref_succeed_ITI_err=[sum(pref_succeed_iti_err.*len),...
                                sum((onechain.trials(find(pref_delay_err_trls)+1,1)-onechain.trials(pref_delay_err_trls,2))./sps-3)]; %rwd + test

                            % succeed ITI nonpref
                            nonpref_succeed_iti=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp_pref ... % WT, nonpref
                                & trl_align(:,2)>=(trl_align(:,4)+4)...  % not in delay or test
                                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14)); % 1s samp, 1s test, 2s rwd
                            freqstats.nonpref_succeed_ITI=[sum(nonpref_succeed_iti.*len),...
                                sum((onechain.trials(find(nonpref_delay_trls)+1,1)-onechain.trials(nonpref_delay_trls,2))./sps-3)];
                        end
                        onechain.trials(end,:)=[];

                        % precede preferred, non preferred

                        % precede ITI pref correct
                        onechain.trials=[repmat(onechain.trials(1,1)-14*sps,1,10);onechain.trials];
                        pref_precede_iti=all(trl_align(:,12:14)==1,2)...
                            & (trl_align(:,2)>=(trl_align(:,4)+4) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
                        freqstats.pref_precede_ITI=[sum(pref_precede_iti.*len),sum((onechain.trials(find(pref_delay_trls)+1,1)-onechain.trials(pref_delay_trls,2))./sps-3)]; %rwd + test
                        if ~nnonmem
                            % precede ITI pref error
                            pref_succeed_iti_err=trl_align(:,13)==0 & trl_align(:,14)==1 ...
                                & (trl_align(:,2)>=(trl_align(:,4)+4) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
                            freqstats.pref_precede_ITI_err=[sum(pref_succeed_iti_err.*len),sum((onechain.trials(find(pref_delay_err_trls)+1,1)-onechain.trials(pref_delay_err_trls,2))./sps-3)]; %rwd + test

                            %  precede ITI nonpref
                            nonpref_precede_iti=all(trl_align(:,12:13)==1,2) & trl_align(:,10)~=samp_pref...
                                & (trl_align(:,2)>=(trl_align(:,4)+4) |(trl_align(:,2)==-1 & trl_align(:,9)<11)); % other trial | first trial
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

                        motif_trl.(dd{1}).(ww{1}).(cc{1}).trl_align=trl_align;
                        motif_trl.(dd{1}).(ww{1}).(cc{1}).freqstats=freqstats;
                    end
                end
            end

            sum_stats=[];
            raw=cell2struct({[];[];cell(0);cell(0)},{'count','time','condition','tag'});
            for dd=reshape(fieldnames(motif_trl),1,[])
                for ww=reshape(fieldnames(motif_trl.(dd{1})),1,[])
                    for cc=reshape(fieldnames(motif_trl.(dd{1}).(ww{1})),1,[])
                        sum_stats=[sum_stats,cellfun(@(x) x(1)./x(2),struct2cell(motif_trl.(dd{1}).(ww{1}).(cc{1}).freqstats))];
                        raw.count=[raw.count,cellfun(@(x) x(1),struct2cell(motif_trl.(dd{1}).(ww{1}).(cc{1}).freqstats))];
                        raw.time=[raw.time,cellfun(@(x) x(2),struct2cell(motif_trl.(dd{1}).(ww{1}).(cc{1}).freqstats))];
                        sess_str=regexp(cc{1},'^s\d{1,3}(?=(c|r))','match','once');
                        raw.condition{end+1}=[sess_str,replace(ww{1},'olf_s','o'),dd{1}];
                        raw.tag{end+1}=regexp(cc{1},'(c|r)\d*$','match','once');
                    end
                end
            end
        end

        %%
        function rec_dur=sessid2length(sessidx)
            persistent sessidx_ rec_dur_
            if isempty(sessidx_) || sessidx~=sessidx_
                metaf=dir(fullfile(...
                    ephys.util.getHomedir('type','raw'),...
                    replace(ephys.sessid2path(sessidx,'criteria','WT'),'\',filesep()),...
                    "*.ap.meta"));
                ts=textscan(...
                    fopen(fullfile(metaf.folder,metaf.name)),...
                    '%s','Delimiter',{'\n'});
                rec_dur=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''))...
                    ./385/2;
                rec_dur_=rec_dur;
                sessidx_=sessidx;
            else
                rec_dur=rec_dur_;
            end
        end

        %%
        function [fhb,fhs]=plot_replay(stats,xlbl,opt)
            arguments
                stats double
                xlbl char
                opt.title (1,:) char = 'chains'
            end

            fhb=figure();
            boxplot(stats.','Colors','k','Symbol','c.')
            ylim([0,1.5])
            set(gca(),'XTick',1:size(stats,1),'XTickLabel',xlbl,'YScale','linear')
            for jj=2:size(stats,1)
                pp=ranksum(stats(1,:),stats(jj,:));
                text(jj,1.5,sprintf('%.3f',pp),'VerticalAlignment','top','HorizontalAlignment','center')
            end
            title(opt.title)
            ylabel('Motif spike frequency (Hz)')

            fhs=figure();
            hold on
            for ii=1:size(stats,1)
                swarmchart(repmat(ii,size(stats,2),1),stats(ii,:),16,'.')
            end
            ylim([-0.1,1.5])
            set(gca(),'XTick',1:size(stats,1),'XTickLabel',xlbl,'YScale','linear','TickLabelInterpreter','none')
            title(opt.title)
            ylabel('Motif spike frequency (Hz)')
        end

        function [per_sess_struct,per_sess_mat]=stats_replay_sess(stats_all,opt)
            arguments
                stats_all cell
                opt.feat_sel = []
            end

            per_sess_struct=struct();
            for ii=1:numel(stats_all)
                motif_set=stats_all{ii};
                if isempty(opt.feat_sel)
                    featsel=true(size(motif_set.count(:,1)));
                else
                    featsel=opt.feat_sel;
                end
                for jj=1:size(motif_set.count,2)
                    onekey=motif_set.condition{jj};
                    if ~isfield(per_sess_struct,onekey)
                        per_sess_struct.(onekey)=cell2struct({motif_set.count(featsel,jj).';motif_set.time(featsel,jj).';motif_set.tag(jj)},{'count';'time';'tag'});
                    else
                        if isequaln(per_sess_struct.(onekey).time,motif_set.time(featsel,jj).')
                            per_sess_struct.(onekey).count=[per_sess_struct.(onekey).count;motif_set.count(featsel,jj).'];
                            per_sess_struct.(onekey).tag=[per_sess_struct.(onekey).tag;motif_set.tag(jj)];
                        else
                            disp("time does not match")
                            keyboard()
                        end
                    end
                end
            end

            outcell=struct2cell(per_sess_struct);
            per_sess_mat=cell2mat(cellfun(@(x) sum(x.count,1)./x.time,outcell,'UniformOutput',false));
        end


        function fhb=plot_replay_sess(per_sess_mat,xlbl,opt)
            arguments
                per_sess_mat double
                xlbl char
                opt.title (1,:) char = 'chains'
                opt.ref_line (1,1) logical = false
                opt.ref_p_value (1,1) logical = true
                opt.median_value (1,1) logical = false
            end
            cmatv=reshape(per_sess_mat,[],1);
            cmatg=reshape(repmat(1:size(per_sess_mat,2),size(per_sess_mat,1),1),[],1);


            fhb=figure('Position',[100,100,600,400]);
            boxplot(cmatv(isfinite(cmatv))+eps,cmatg(isfinite(cmatv)),'Colors','k','Symbol','c.')
            if opt.ref_line
                yline(median(per_sess_mat(:,1)),'--r');
            end
            set(gca(),'XTick',1:size(per_sess_mat,2),'XTickLabel',xlbl,'XTickLabelRotation',90,'YScale','log')
            if opt.ref_p_value
                for jj=2:size(per_sess_mat,2)
                    finisel=isfinite(per_sess_mat(:,jj));
                    pp=signrank(per_sess_mat(finisel,1),per_sess_mat(finisel,jj));
                    text(jj,0.01,sprintf('%.3f',pp),'VerticalAlignment','bottom','HorizontalAlignment','center')
                end
            end

            if opt.median_value
                for jj=1:size(per_sess_mat,2)
                    mm=median(cmatv(cmatg==jj & isfinite(cmatv)));
                    text(jj,1,sprintf('%.1f',mm),'VerticalAlignment','bottom','HorizontalAlignment','center')
                end
            end

            ylim([0.002,300]);
            %
            % fhs=figure();
            % hold on
            % for ii=1:size(stats_all,1)
            %     swarmchart(repmat(ii,size(stats_all,2),1),stats_all(ii,:),16,'.')
            % end
            % ylim([-0.1,1.5])
            % set(gca(),'XTick',1:size(stats_all,1),'XTickLabel',xlbl,'YScale','linear','TickLabelInterpreter','none')
            title(opt.title)
            ylabel('Motif spike frequency (Hz)')
        end

        function fhb=plot_replay_sess_ci(per_sess_mat,xlbl,opt)
            arguments
                per_sess_mat double
                xlbl char
                opt.title (1,:) char = 'chains'
                opt.ref_line (1,1) logical = false
                opt.ref_p_value (1,1) logical = true
                opt.median_value (1,1) logical = false
            end
            cmatv=reshape(per_sess_mat,[],1);
            cmatg=reshape(repmat(1:size(per_sess_mat,2),size(per_sess_mat,1),1),[],1);

            mm=nanmedian(per_sess_mat);
            ci=bootci(1000,@(x) nanmedian(x),per_sess_mat);

            fhb=figure('Position',[100,100,600,400]);
            hold on
            bar(mm.','grouped','FaceColor','none','EdgeColor','k')
            errorbar(1:numel(mm),mm,ci(1,:)-mm,ci(2,:)-mm,'k.');
            if opt.ref_line
                yline(median(per_sess_mat(:,1)),'--r');
            end
            set(gca(),'XTick',1:size(per_sess_mat,2),'XTickLabel',xlbl,'XTickLabelRotation',90,'YScale','log')
            if opt.ref_p_value
                for jj=2:size(per_sess_mat,2)
                    finisel=isfinite(per_sess_mat(:,jj));
                    pp=signrank(per_sess_mat(finisel,1),per_sess_mat(finisel,jj));
                    text(jj,0.01,sprintf('%.3f',pp),'VerticalAlignment','bottom','HorizontalAlignment','center')
                end
            end

            if opt.median_value
                for jj=1:size(per_sess_mat,2)
                    mm=median(cmatv(cmatg==jj & isfinite(cmatv)));
                    text(jj,1,sprintf('%.1f',mm),'VerticalAlignment','bottom','HorizontalAlignment','center')
                end
            end

            ylim([0.1,10]);
            title(opt.title)
            ylabel('Motif spike frequency (Hz)')
        end


        function fhb=plot_replay_cross_sess(cross_sess_mat,xlbl,opt)
            arguments
                cross_sess_mat cell
                xlbl char
                opt.title (1,:) char = 'chains'
                opt.median_value (1,1) logical = false
            end
            mergemat=cell2mat(arrayfun(@(x) [cross_sess_mat{x},repmat(x,numel(cross_sess_mat{x}),1)],(1:numel(cross_sess_mat)).','UniformOutput',false));
            mergemat=mergemat(isfinite(mergemat(:,1)),:);

            fhb=figure('Position',[100,100,600,400]);
            for ii=1:size(cross_sess_mat,2)
                boxplot(mergemat(:,1),mergemat(:,2),'Colors','k','Symbol','c.')
            end
            set(gca(),'XTick',1:size(cross_sess_mat,2),'XTickLabel',xlbl,'XTickLabelRotation',90,'YScale','log')

            if opt.median_value
                for jj=1:size(cross_sess_mat,2)
                    mm=median(cross_sess_mat{jj}(isfinite(cross_sess_mat{jj})));
                    text(jj,1,sprintf('%.1f',mm),'VerticalAlignment','bottom','HorizontalAlignment','center')
                end
            end

            ylim([0.002,300]);
            %
            % fhs=figure();
            % hold on
            % for ii=1:size(stats_all,1)
            %     swarmchart(repmat(ii,size(stats_all,2),1),stats_all(ii,:),16,'.')
            % end
            % ylim([-0.1,1.5])
            % set(gca(),'XTick',1:size(stats_all,1),'XTickLabel',xlbl,'YScale','linear','TickLabelInterpreter','none')
            title(opt.title)
            ylabel('Motif spike frequency (Hz)')
        end


        function showcase()
            %%
            global_init;
            load(fullfile(gather_config.odpath,'Tempdata','TEMP230602.mat'),'sschain_trl');
            [chain_replay,~,~]=wave.replay.stats(sschain_trl,'var_len',false);
            %%
            for sampkey=["olf_s1","olf_s2"]
                fns=intersect(fieldnames(chain_replay.d3.(sampkey)),fieldnames(chain_replay.d6.(sampkey)));
                durkey="d3";%for durkey=["d3","d6"]
                for fnidx=1:numel(fns)
                    onechain=chain_replay.(durkey).(sampkey).(fns{fnidx});
                    onechain.trials=[onechain.trials,(1:size(onechain.trials,1)).'];
                    onechain.trl_align=unique(onechain.trl_align,'rows');
                    if ~any(onechain.trl_align(:,7),'all')
                        continue
                    end
                    pref_samp=onechain.trl_align(find(onechain.trl_align(:,7)==1,1),3);
                    % ===================
                    trl_pool=onechain.trials(onechain.trials(:,9)==1,:);
                    for tidx=1:(size(trl_pool,1)-15)
                        if numel(unique(trl_pool(tidx:(tidx+3),8)))>1
                            continue
                        end
                        % survey for delay selectivity
                        if nnz(ismember(...
                                onechain.trl_align(...
                                onechain.trl_align(:,2)<(1+onechain.trl_align(:,4))... % in delay
                                & onechain.trl_align(:,7)==1 ... % prefered
                                ,1),... % extract trialnum
                                trl_pool(tidx:(tidx+15),11)...
                                ))<10
                            continue
                        end
                        % exclude nonspecific activity
                        if nnz(ismember(...
                                onechain.trl_align(...
                                onechain.trl_align(:,2)<(1+onechain.trl_align(:,4))... % in delay
                                & onechain.trl_align(:,7)==0 ... % prefered
                                ,1),... % extract trialnum
                                trl_pool(tidx:(tidx+15),11)...
                                ))>4
                            continue
                        end

                        wt_trls=trl_pool(tidx:(tidx+15),:);
                        sps=30000;
                        sbound=[onechain.trials(1,1),onechain.trials(end,2)];
                        hexcmap=["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];

                        sessid=str2double(regexp(fns{fnidx},'(?<=^s)\d{1,3}(?=(c|r))','match','once'));
                        [SPKID,SPKTS,~,~,~,~]=ephys.getSPKID_TS(sessid,'keep_trial',false);

                        % ===================

                        fh=figure('Position',[1280,320,1280,480]);
                        tiledlayout(3,1)
                        % before session >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        nexttile(1)
                        hold on
                        

                        before_sess_sel=onechain.ts(:,end)<sbound(1) & onechain.ts(:,1)>(sbound(1)-660*sps);
                        bg_fr=nan(size(onechain.ts,2),1);
                        for yy=1:size(onechain.ts,2)
                            suid=onechain.meta{1}(yy);
                            suspk=SPKTS(SPKID==suid & SPKTS<sbound(1) & SPKTS>(sbound(1)-660*sps));
                            bg_fr(yy)=numel(suspk)./660;
                            plot((onechain.ts(before_sess_sel,yy)-sbound(1))./sps,yy,'|','Color',hexcmap(yy));
                        end
                        xlim([-300,0])
                        ylim([0.5,yy+0.5])
                        set(gca(),YDir="reverse")
                        % set(gca(),'YTick',10:5:29,'YTickLabel',100:50:290,'XTick',0:5:10,'XTickLabel',10:-5:0)
                        xlabel('Time relative to session onset (s)')

                        % during session >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        nexttile(2)
                        hold on
                        %
                        % sample, delay n' test
                        yspan=size(onechain.ts,2)+0.5;
                        for ii=1:size(wt_trls,1)
                            fill(([wt_trls(ii,1),wt_trls(ii,1)+sps,wt_trls(ii,1)+sps,wt_trls(ii,1)]-sbound(1))./sps,[0.5,0.5,yspan,yspan],'k','EdgeColor','none','FaceAlpha',0.1);
                            fill(([wt_trls(ii,2),wt_trls(ii,2)+sps,wt_trls(ii,2)+sps,wt_trls(ii,2)]-sbound(1))./sps,[0.5,0.5,yspan,yspan],'b','EdgeColor','none','FaceAlpha',0.1);
                            if wt_trls(ii,5)==pref_samp
                                fill(([wt_trls(ii,1)+sps,wt_trls(ii,2),wt_trls(ii,2),wt_trls(ii,1)+sps]-sbound(1))./sps,[0.5,0.5,yspan,yspan],'r','EdgeColor','none','FaceAlpha',0.1);
                            else
                                fill(([wt_trls(ii,1)+sps,wt_trls(ii,2),wt_trls(ii,2),wt_trls(ii,1)+sps]-sbound(1))./sps,[0.5,0.5,yspan,yspan],'y','EdgeColor','none','FaceAlpha',0.1);
                            end
                        end
                        % motif onsets
                        bound=[wt_trls(1,1),wt_trls(end,2)]+sps.*[-5,10];
                        motif_spks=onechain.ts((onechain.ts(:,1)>=bound & onechain.ts(:,end)<bound(2)),:);

                        for yy=1:size(motif_spks,2)
                            plot((motif_spks(:,yy)-sbound(1))./sps,yy,'|','Color',hexcmap(yy));
                        end
                        xlim((bound-sbound(1))./sps)
                        xlabel("Time since session onset (s)")
                        ylim([0.5,yy+0.5])
                        set(gca(),YDir="reverse")
                        %
                        % after session >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        nexttile(3)
                        hold on
                        after_sess_sel=onechain.ts(:,1)>sbound(2) & onechain.ts(:,1)<(sbound(2)+660*sps);
                        for yy=1:size(onechain.ts,2)
                            plot((onechain.ts(after_sess_sel,yy)-sbound(2))./sps,yy,'|','Color',hexcmap(yy));
                        end
                        
                        ylim([0.5,yy+0.5])
                        set(gca(),YDir="reverse")
                        % set(gca(),'YTick',10:5:29,'YTickLabel',100:50:290,'XTick',0:5:10,'XTickLabel',10:-5:0)
                        xlabel('Time since last trial (s)')
                        xlim([60,360])
                        % keyboard()
                        sgtitle(durkey+","+sampkey+","+fns{fnidx}+","+tidx,Interpreter='none')
                        disp(bg_fr);
                    end
                end
            end
            %% =============================================
            su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
            for ii=sschain_trl.d3.olf_s2.s102c1648.meta{1}
                disp(su_meta.reg_tree(5,su_meta.sess==102 & su_meta.allcid==ii))
            end


        end

        function [fh,out]=region_replay(motif_replay,opt)
            arguments
                motif_replay
                opt.reg="HIP"
            end
            out=cell2struct({[];[]},{'REG','Others'});
            su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
            for dd=reshape(string(fieldnames(motif_replay)),1,[])
                for samp=reshape(string(fieldnames(motif_replay.(dd))),1,[])
                    for condi=reshape(string(fieldnames(motif_replay.(dd).(samp))),1,[])
                        % sessid=str2double(regexp(condi,'(?<=^s)\d{1,3}(?=(c|r))','match','once'));
                        onemotif=motif_replay.(dd).(samp).(condi);
                        sessid=onemotif.meta{1};
                        cids=onemotif.meta{2};
                        motif_reg=arrayfun(@(x) string(su_meta.reg_tree{5,su_meta.sess==sessid & su_meta.allcid==x}),cids);
                        onefreq=[onemotif.freqstats.pref_delay_correct,...
                                onemotif.freqstats.pref_succeed_ITI,...
                                onemotif.freqstats.before_session,...
                                onemotif.freqstats.after_session];

                        if any(motif_reg==opt.reg)
                            out.REG=[out.REG;onefreq(1:2:7)./onefreq(2:2:8)];
                        else
                            out.Others=[out.Others;onefreq(1:2:7)./onefreq(2:2:8)];
                        end
                    end
                end
            end
            dd=[out.REG;out.Others];
            gg=[repmat([1:2:8],size(out.REG,1),1);...
                repmat([2:2:8],size(out.Others,1),1)];
            fh=figure();
            boxplot(reshape(dd,[],1),reshape(gg,[],1),'Symbol','c.','Colors','k')
            set(gca,'XTick',2:2:8,'XTickLabel',{'Delay','aft.ITI','Before sess.','After Sess.'})
            ylim([0,2])

            for ii=1:2:8
                pp=ranksum(dd(gg==ii), dd(gg==ii+1));
                text(ii+0.5,2,sprintf('%.3f',pp));
            end

        end



        function exception_check()
            [spkID,spkTS,trials,SU_id,folder,FT_SPIKE]=ephys.getSPKID_TS(106,'keep_trial',true);
            errsel=trials(:,10)==0;
            sampsel=trials(:,5)==8;
            dur3sel=trials(:,8)==3;
            dur6sel=trials(:,8)==6;
            susel=strcmp(FT_SPIKE.label,'10202');

            corrct3=nnz(... % trial
                ismember(FT_SPIKE.trial{susel},find(~errsel & sampsel & dur3sel))...
                ... % time
                & FT_SPIKE.time{susel}>=1 & FT_SPIKE.time{susel}<4)...
                ... % total delay time
                ./(nnz(~errsel & sampsel & dur3sel).*3)

            error3=nnz(... % trial
                ismember(FT_SPIKE.trial{susel},find(errsel & sampsel & dur3sel))...
                ... % time
                & FT_SPIKE.time{susel}>=1 & FT_SPIKE.time{susel}<4)...
                ... % total delay time
                ./(nnz(errsel & sampsel & dur3sel).*3)


            corrct6=nnz(... % trial
                ismember(FT_SPIKE.trial{susel},find(~errsel & sampsel & dur3sel))...
                ... % time
                & FT_SPIKE.time{susel}>=1 & FT_SPIKE.time{susel}<7)...
                ... % total delay time
                ./(nnz(~errsel & sampsel & dur6sel).*6)

            error6=nnz(... % trial
                ismember(FT_SPIKE.trial{susel},find(errsel & sampsel & dur3sel))...
                ... % time
                & FT_SPIKE.time{susel}>=1 & FT_SPIKE.time{susel}<7)...
                ... % total delay time
                ./(nnz(errsel & sampsel & dur6sel).*6)


            for ii=1:size(trials,1)
                if rem(ii,40)==1
                    figure()
                    hold on
                end
                if all(trials(ii,9:10)==1,2)
                    plot(FT_SPIKE.time{susel}(FT_SPIKE.trial{susel}==ii),ii,'k|')
                elseif trials(ii,10)==0
                    plot(FT_SPIKE.time{susel}(FT_SPIKE.trial{susel}==ii),ii,'r|')
                else
                    plot(FT_SPIKE.time{susel}(FT_SPIKE.trial{susel}==ii),ii,'|','Color',[0.5,0.5,0.5])
                end
            end
        end



    end
end