classdef delay_vs_iti_per_sec < handle
    methods (Static)
        % function [cover_per_sec,fh]=run(covered,trials_dict,opt)
        %     arguments
        %         covered = []
        %         trials_dict = []
        %         opt.skip_save (1,1) logical = false
        %     end
        %     if isempty(covered)
        %         fstr=load(fullfile('binary','delay_iti_runlength_covered.mat'),'covered_tbl');
        %         covered=fstr.covered_tbl;
        %         clear fstr;
        %     end
        %     if isempty(trials_dict)
        %         load(fullfile('binary','trials_dict.mat'),'trials_dict');
        %     end
        % 
        %     [cover_per_sec,mdm,ci,per_sess]=wave.replay.delay_vs_iti_per_sec.statsOne(covered,trials_dict);
        %     fh=wave.replay.delay_vs_iti_per_sec.plotOne(mdm,ci,per_sess);
        %     if ~opt.skip_save
        %         blame=vcs.blame();
        %         save(fullfile('binary','motif_cover_per_sec.mat'),'cover_per_sec','blame');
        %     end
        % end

        function [cover_per_sec,fh]=run_shuf(covered,trials_dict,opt)
            arguments
                covered = []
                trials_dict = []
                opt.skip_save (1,1) logical = false
            end
            if isempty(trials_dict)
                load(fullfile('binary','trials_dict.mat'),'trials_dict');
            end
            if isempty(covered)
                fstr=load(fullfile('binary','delay_iti_runlength_covered.mat'),'covered_tbl');
            end
            cover_per_sec=wave.replay.delay_vs_iti_per_sec.statsOne(fstr.covered_tbl,trials_dict);
            if ~opt.skip_save
                blame=vcs.blame();
                save(fullfile('binary','motif_cover_per_sec.mat'),'cover_per_sec','blame');
            end
            shufstr=load(fullfile("binary","delay_iti_runlength_covered_shuf.mat"),'covered_shuf');
            fh=wave.replay.delay_vs_iti_per_sec.plotshuf(cover_per_sec,shufstr.covered_shuf);
        end

        function [cover_per_sec,mdm,ci,per_sess]=statsOne(covered,trials_dict)
            sps=30000;
            % per session
            cover_per_sec=[];
            for sessidx=1:size(covered,1)
                sess=covered.session(sessidx);
                switch covered.wave(sessidx)
                    case "s1d3"
                        samp=4;delay=3;
                    case "s1d6"
                        samp=4;delay=6;
                    case "s2d3"
                        samp=8;delay=3;
                    case "s2d6"
                        samp=8;delay=6;
                end

                disp(sess)
                session_tick=wave.replay.sessid2length(sess);
                trials=cell2mat(trials_dict(sess));
                surround_sec=(trials(1,1)./sps-60)+((session_tick-trials(end,2))./sps-60-1);
                surround_motif_sec=nnz(covered.out_task{sessidx})./10000;

                %  corresponding network in pre task, post task

                % per preferred trial
                trl_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
                np_trl_sel=find(trials(:,5)==setdiff([4,8],samp) & trials(:,8)==delay & all(trials(:,9:10)>0,2));

                pref_delay_sec=sum(diff(trials(trl_sel,1:2),1,2)./sps-1);
                delayc=covered.delay{sessidx};
                delay_motif_sec=0;
                for tt=reshape(trl_sel,1,[])
                    delay_motif_sec=delay_motif_sec+nnz(delayc(floor(trials(tt,1)./3):ceil(trials(tt,2)./3)))./10000;
                end

                np_delay_sec=sum(diff(trials(np_trl_sel,1:2),1,2)./sps-1);
                npdelayc=covered.npdelay{sessidx};
                npdelay_motif_sec=0;
                for tt=reshape(np_trl_sel,1,[])
                    npdelay_motif_sec=npdelay_motif_sec+nnz(npdelayc(floor(trials(tt,1)./3):ceil(trials(tt,2)./3)))./10000;
                end

                trials(end+1,:)=min(session_tick-3,trials(end,2)+14*sps);
                pref_succeed_iti_sec=sum((trials(trl_sel+1,1)-trials(trl_sel,2))./sps-4); % 1s test + 3s response
                itic=covered.iti{sessidx};
                iti_motif_sec=0;
                for tt=reshape(trl_sel,1,[])
                    iti_motif_sec=iti_motif_sec+nnz(itic(floor(trials(tt,2)./3):ceil(trials(tt+1,1)./3)))./10000;
                end

                npiti_sec=sum((trials(np_trl_sel+1,1)-trials(np_trl_sel,2))./sps-4); % 1s test + 3s response
                npitic=covered.npiti{sessidx};
                npiti_motif_sec=0;
                for tt=reshape(np_trl_sel,1,[])
                    npiti_motif_sec=npiti_motif_sec+nnz(npitic(floor(trials(tt,2)./3):ceil(trials(tt+1,1)./3)))./10000;
                end

                trials(end,:)=[];

                cover_per_sec=[cover_per_sec;...
                    cell2table({sess,samp,delay,[delay_motif_sec,pref_delay_sec],[iti_motif_sec,pref_succeed_iti_sec],[surround_motif_sec,surround_sec],[npdelay_motif_sec,np_delay_sec],[npiti_motif_sec,npiti_sec]},...
                    'VariableNames',{'session','sample','dur','delay','iti','surround','npdelay','npiti'})];
            end
            dps=cover_per_sec.delay;
            dprop=(dps(:,1)./dps(:,2));

            npdps=cover_per_sec.npdelay;
            npdprop=(npdps(:,1)./npdps(:,2));

            ips=cover_per_sec.iti;
            iprop=(ips(:,1)./ips(:,2));

            npips=cover_per_sec.npiti;
            npiprop=(npips(:,1)./npips(:,2));

            sps=cover_per_sec.surround;
            sprop=(sps(:,1)./sps(:,2)); % Due to 4 waves

            mdm=median([dprop,npdprop,iprop,npiprop,sprop]);
            ci=bootci(1000,@(x) median(x),[dprop,npdprop,iprop,npiprop,sprop]);
            per_sess=cell2struct({dprop;npdprop;iprop;npiprop;sprop},{'delay','iti','npdelay','npiti','offtask'});
        end
        function fh=plotOne(mdm,ci,per_sess)

            %% plot
            % mm=[mean(dprop),mean(npdprop),mean(iprop),mean(npiprop),mean(sprop),];
            % sem=[std(dprop),std(npdprop),std(iprop),std(npiprop),std(sprop)]./sqrt([numel(dprop),numel(npdprop),numel(iprop),numel(npiprop),numel(sprop)]);

            fh=figure('Position',[100,100,400,300]);
            hold on
            bh=bar(mdm.*100,'FaceColor','none','FaceColor','k');
            errorbar(bh.XEndPoints,bh.YEndPoints,(ci(1,:)-mdm).*100,(ci(2,:)-mdm).*100,'k.')
            set(gca,'XTick',1:5,'XTickLabelRotation',90,'XTickLabel',{'Delay','NPDelay','ITI','NPITI','Surround'});
            ylabel('Motif duration / total duration (%)');

            pnpdelay=signrank(per_sess.delay,per_sess.npdelay);
            piti=signrank(per_sess.delay,per_sess.iti);
            poff=signrank(per_sess.delay,per_sess.offtask);
            piti_iti=signrank(per_sess.iti,per_sess.npiti);

            title(sprintf('%s%.4f,','p-np',pnpdelay,'d-iti',piti,'iti-iti',piti_iti,'otask',poff));
            savefig(fh,fullfile('binary','delay_vs_iti_per_sec.fig'));
            appendfig('close',true,'tag','delay vs iti composite motif coverage time, delay_vs_iti_per_sec.m');
        end

        function fh=plotshuf(observ,shuf)
            dps=observ.delay;
            dprop=(dps(:,1)./dps(:,2));
            npdps=observ.npdelay;
            npdprop=(npdps(:,1)./npdps(:,2));
            ips=observ.iti;
            iprop=(ips(:,1)./ips(:,2));
            sps=observ.surround;
            soff=(sps(:,1)./sps(:,2)); % Due to 4 waves
            mdm=median([dprop,npdprop,iprop,soff]);
            ci=bootci(1000,@(x) median(x),[dprop,npdprop,iprop,soff]);
            %% shuffle
            shufdps=shuf.delay;
            shufdprop=(shufdps(:,1)./shufdps(:,2));
            shufnpdps=shuf.npdelay;
            shufnpdprop=(shufnpdps(:,1)./shufnpdps(:,2));
            shufips=shuf.iti;
            shufiprop=(shufips(:,1)./shufips(:,2));
            shufoffps=shuf.surround;
            shufoff=(shufoffps(:,1)./shufoffps(:,2)); % Due to 4 waves
            shufmdm=median([shufdprop,shufnpdprop,shufiprop,shufoff]);
            shufci=bootci(1000,@(x) median(x),[shufdprop,shufnpdprop,shufiprop,shufoff]);

            fh=figure('Position',[100,100,400,300]);
            hold on
            bh=bar(100.*[mdm;shufmdm].','grouped','FaceColor','k');
            bh(2).FaceColor='w';
            errorbar(bh(1).XEndPoints,bh(1).YEndPoints,(ci(1,:)-mdm).*100,(ci(2,:)-mdm).*100,'k.')
            errorbar(bh(2).XEndPoints,bh(2).YEndPoints,(shufci(1,:)-shufmdm).*100,(shufci(2,:)-shufmdm).*100,'k.')
            set(gca,'XTick',1:5,'XTickLabelRotation',90,'XTickLabel',{'Delay','NPDelay','ITI','Off-task'});
            ylabel('Motif duration / total duration (%)');
            
            pnpdelay=signrank(dprop,npdprop);
            piti=signrank(dprop,iprop);
            psur=signrank(dprop,soff);
            
            zscores=(mean([dprop,npdprop,iprop,soff])-mean([shufdprop,shufnpdprop,shufiprop,shufoff]))./std([shufdprop,shufnpdprop,shufiprop,shufoff]);
            title({sprintf('%s%.4f,','p-np',pnpdelay,'d-iti',piti,'otask',psur);...
                sprintf('z%.3f',zscores)});
            savefig(fh,fullfile('binary','delay_vs_iti_per_sec_w_shuf.fig'));
            appendfig('close',isunix,'tag','delay vs iti composite motif coverage time with shuffle, delay_vs_iti_per_sec.m');        end

    end
end
