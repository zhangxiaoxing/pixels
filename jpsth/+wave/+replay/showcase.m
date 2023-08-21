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

