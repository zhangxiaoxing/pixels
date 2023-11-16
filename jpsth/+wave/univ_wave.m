classdef univ_wave < handle
    % TODO: universal wave
    methods (Static)
        function plot_w_delay(opt)
            arguments
                opt.itisort (1,1) logical = false
            end
            delay=6;
            global_init;
            wrs_mux_meta=ephys.get_wrs_mux_meta();
            com_map=wave.get_pct_com_map(wrs_mux_meta,'odor_only',true);
            [~,imdata]=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0.1,0.7],'gauss2d',true,'delay',delay);
            [imdata.s1n.s1iti,imdata.s1n.s2iti,imdata.s2n.s1iti,imdata.s2n.s2iti]=deal([]);

            sessid=-1;
            for ii=1:numel(imdata.s1n.sess)
                if ~(sessid==imdata.s1n.sess(ii))
                    sessid=imdata.s1n.sess(ii);
                    warning(string(sessid))
                    [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
                    trls=FT_SPIKE.trialinfo;
                    s1trl=find(trls(:,5)==4 & all(trls(:,9:10),2) & trls(:,8)==delay);
                    s2trl=find(trls(:,5)==8 & all(trls(:,9:10),2) & trls(:,8)==delay);
                end

                cid=imdata.s1n.id(ii);
                susel=strcmp(FT_SPIKE.label,num2str(cid));
                s1tsel=ismember(FT_SPIKE.trial{susel},s1trl)...
                    & FT_SPIKE.time{susel} >= -6+delay...
                    & FT_SPIKE.time{susel} < 10+delay;

                s1tfreq=histcounts(FT_SPIKE.time{susel}(s1tsel)+6-delay,0:0.25:16)./numel(s1trl);
                imdata.s1n.s1iti=[imdata.s1n.s1iti;s1tfreq];

                s2tsel=ismember(FT_SPIKE.trial{susel},s2trl)...
                    & FT_SPIKE.time{susel} >= -6+delay...
                    & FT_SPIKE.time{susel} < 10+delay;
                s2tfreq=histcounts(FT_SPIKE.time{susel}(s2tsel)+6-delay,0:0.25:16)./numel(s2trl);
                imdata.s1n.s2iti=[imdata.s1n.s2iti;s2tfreq];

            end

            sessid=-1;
            for ii=1:numel(imdata.s2n.sess)
                if ~(sessid==imdata.s2n.sess(ii))
                    sessid=imdata.s2n.sess(ii);
                    warning(string(sessid))
                    [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
                    trls=FT_SPIKE.trialinfo;
                    s1trl=find(trls(:,5)==4 & all(trls(:,9:10),2) & trls(:,8)==delay);
                    s2trl=find(trls(:,5)==8 & all(trls(:,9:10),2) & trls(:,8)==delay);
                end

                cid=imdata.s2n.id(ii);
                susel=strcmp(FT_SPIKE.label,num2str(cid));
                s1tsel=ismember(FT_SPIKE.trial{susel},s1trl)...
                    & FT_SPIKE.time{susel} >= -6+delay...
                    & FT_SPIKE.time{susel} < 10+delay;
                s1tfreq=histcounts(FT_SPIKE.time{susel}(s1tsel)+6-delay,0:0.25:16)./numel(s1trl);
                imdata.s2n.s1iti=[imdata.s2n.s1iti;s1tfreq];

                s2tsel=ismember(FT_SPIKE.trial{susel},s2trl)...
                    & FT_SPIKE.time{susel} >= -6+delay...
                    & FT_SPIKE.time{susel} < 10+delay;
                s2tfreq=histcounts(FT_SPIKE.time{susel}(s2tsel)+6-delay,0:0.25:16)./numel(s2trl);
                imdata.s2n.s2iti=[imdata.s2n.s2iti;s2tfreq];

            end


            S=max(imdata.s1n.s1iti(:,5:28)-s1itimm,[],2);
            s1s1t=(imdata.s1n.s1iti-s1itimm)./S;
            s1s2t=(imdata.s1n.s2iti-s1itimm)./S;

            s2itimm=mean((imdata.s2n.s1iti(:,5:28)+imdata.s2n.s2iti(:,5:28))./2,2);
            S=max(imdata.s2n.s2iti(:,5:28)-s2itimm,[],2);
            s2s1t=(imdata.s2n.s1iti-s2itimm)./S;
            s2s2t=(imdata.s2n.s2iti-s2itimm)./S;

            if opt.itisort
                % alternative sorts
                % ##########################################################
                s1s1tc=s1s1t(:,45:64);
                s1s1tc(s1s1tc<0)=0;
                localcoms1=sum(repmat(1:20,size(s1s1t,1),1).*s1s1tc,2)./sum(s1s1tc,2); % 1+6+1+3
                localcoms1(all(s1s1tc==0,2))=-1;
                [~,localsorts1]=sort(localcoms1);

                s2s2tc=s2s2t(:,45:64);
                s2s2tc(s2s2tc<0)=0;
                localcoms2=sum(repmat(1:20,size(s2s2t,1),1).*s2s2tc,2)./sum(s2s2tc,2); % 1+6+1+3
                localcoms2(all(s2s2tc==0,2))=-1;
                [~,localsorts2]=sort(localcoms2);
                immat={s1s1t(localsorts1,:),s1s2t(localsorts1,:),...
                    s2s1t(localsorts2,:),s2s2t(localsorts2,:)};
            else
                immat={s1s1t(imdata.s1n.sorts,:),s1s2t(imdata.s1n.sorts,:),...
                    s2s1t(imdata.s2n.sorts,:),s2s2t(imdata.s2n.sorts,:)};
            end

            %% PLOT
            scale=[0,0.7];
            gk = fspecial('gaussian', [9 3], 1);
            figure()
            tiledlayout(2,2)
            sidx=1;
            for datum=immat
                nexttile()
                imagesc(conv2(datum{1},gk,'same'),scale);
                set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
                xlim([0.5,size(datum{1},2)+0.5])
                ylim([0.5,size(datum{1},1)+0.5])
                set(gca(),'YDir','normal')
            end
            colormap("turbo")
        end

        function itidata=stats(opt)
            arguments
                opt.rpt (1,1) double =100
                opt.loadfile (1,1) logical = false
                opt.filename (1,:) char = 'temp0808.mat'
            end
            if opt.loadfile
                load(fullfile('binary/',opt.filename),'itidata');
                return
            end
            delay=6;
            global_init;
            wrs_mux_meta=ephys.get_wrs_mux_meta();
            com_map=wave.get_pct_com_map(wrs_mux_meta,'odor_only',true);
            [fh6,imdata]=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0.1,0.7],'gauss2d',true,'delay',delay);
            
            close(fh6.wave4);

            [imdata.s1n.s1iti,imdata.s1n.s2iti,imdata.s2n.s1iti,imdata.s2n.s2iti]=deal([]);

            itidata.s1n.h1h=nan(numel(imdata.s1n.id),100);
            itidata.s1n.h2h=nan(numel(imdata.s1n.id),100);
            itidata.s2n.h1h=nan(numel(imdata.s2n.id),100);
            itidata.s2n.h1h=nan(numel(imdata.s2n.id),100);
            [itidata.s1n.sess,itidata.s1n.cid,itidata.s2n.sess,itidata.s2n.cid]=deal([]);


            sessid=-1;
            for ii=1:numel(imdata.s1n.sess)
                if ~(sessid==imdata.s1n.sess(ii))
                    save(fullfile('binary',opt.filename),'itidata');
                    sessid=imdata.s1n.sess(ii);
                    warning(string(sessid))
                    [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
                    trls=FT_SPIKE.trialinfo;
                    % half-half COM TRLS
                    % ______________________________________________
                    s1trl=find(trls(:,5)==4 & all(trls(:,9:10),2) & trls(:,8)==delay);
                    s2trl=find(trls(:,5)==8 & all(trls(:,9:10),2) & trls(:,8)==delay);
                end

                cid=imdata.s1n.id(ii);
                itidata.s1n.sess=[itidata.s1n.sess;sessid];
                itidata.s1n.cid=[itidata.s1n.cid;cid];
                
                susel=strcmp(FT_SPIKE.label,num2str(cid));
    
                % repeats
                
                for rpt=1:opt.rpt
                    s1trl1h=randsample(s1trl,round(numel(s1trl)./2));
                    s1trl2h=setdiff(s1trl,s1trl1h);

                    s2trl1h=randsample(s2trl,round(numel(s2trl)./2));
                    s2trl2h=setdiff(s2trl,s2trl1h);

                    % one-su com@1h
                    s1tsel1h=ismember(FT_SPIKE.trial{susel},s1trl1h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;

                    s1tfreq1h=histcounts(FT_SPIKE.time{susel}(s1tsel1h)-5-delay,0:0.25:5)./numel(s1trl1h);

                    s2tsel1h=ismember(FT_SPIKE.trial{susel},s2trl1h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;
                    s2tfreq1h=histcounts(FT_SPIKE.time{susel}(s2tsel1h)-5-delay,0:0.25:5)./numel(s2trl1h);

                    s1itimm1h=mean((s1tfreq1h+s2tfreq1h)./2,2);
                    s1s1tc1h=s1tfreq1h-s1itimm1h;
                    s1s1tc1h(s1s1tc1h<0)=0;
                    localcoms1h=sum((1:20).*s1s1tc1h)./sum(s1s1tc1h,2);
                    itidata.s1n.h1h(ii,rpt)=localcoms1h;

                    % TODO one-su com@2h
                    s1tsel2h=ismember(FT_SPIKE.trial{susel},s1trl2h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;

                    s1tfreq2h=histcounts(FT_SPIKE.time{susel}(s1tsel2h)-5-delay,0:0.25:5)./numel(s1trl2h);

                    s2tsel2h=ismember(FT_SPIKE.trial{susel},s2trl2h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;
                    s2tfreq2h=histcounts(FT_SPIKE.time{susel}(s2tsel2h)-5-delay,0:0.25:5)./numel(s2trl2h);

                    s1itimm2h=mean((s1tfreq2h+s2tfreq2h)./2,2);
                    s1s1tc2h=s1tfreq2h-s1itimm2h;
                    s1s1tc2h(s1s1tc2h<0)=0;
                    localcoms2h=sum((1:20).*s1s1tc2h)./sum(s1s1tc2h,2);
                    itidata.s1n.h2h(ii,rpt)=localcoms2h;
                    
                    % if ~all(isfinite([localcoms1h,localcoms2h]),'all')
                    %     keyboard();
                    % end
                end
            end

            sessid=-1;
            for ii=1:numel(imdata.s2n.sess)
                if ~(sessid==imdata.s2n.sess(ii))
                    save(fullfile('binary',opt.filename),'itidata');
                    sessid=imdata.s2n.sess(ii);
                    warning(string(sessid))
                    [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
                    trls=FT_SPIKE.trialinfo;
                    % half-half COM TRLS
                    % ______________________________________________
                    s1trl=find(trls(:,5)==4 & all(trls(:,9:10),2) & trls(:,8)==delay);
                    s2trl=find(trls(:,5)==8 & all(trls(:,9:10),2) & trls(:,8)==delay);
                end

                cid=imdata.s2n.id(ii);
                itidata.s2n.sess=[itidata.s2n.sess;sessid];
                itidata.s2n.cid=[itidata.s2n.cid;cid];

                susel=strcmp(FT_SPIKE.label,num2str(cid));
                % repeats
                
                for rpt=1:opt.rpt
                    s1trl1h=randsample(s1trl,round(numel(s1trl)./2));
                    s1trl2h=setdiff(s1trl,s1trl1h);

                    s2trl1h=randsample(s2trl,round(numel(s2trl)./2));
                    s2trl2h=setdiff(s2trl,s2trl1h);

                    % one-su com@1h
                    s1tsel1h=ismember(FT_SPIKE.trial{susel},s1trl1h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;

                    s1tfreq1h=histcounts(FT_SPIKE.time{susel}(s1tsel1h)-5-delay,0:0.25:5)./numel(s1trl1h);

                    s2tsel1h=ismember(FT_SPIKE.trial{susel},s2trl1h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;
                    s2tfreq1h=histcounts(FT_SPIKE.time{susel}(s2tsel1h)-5-delay,0:0.25:5)./numel(s2trl1h);

                    s1itimm1h=mean((s1tfreq1h+s2tfreq1h)./2,2);
                    s1s1tc1h=s1tfreq1h-s1itimm1h;
                    s1s1tc1h(s1s1tc1h<0)=0;
                    localcoms1h=sum((1:20).*s1s1tc1h)./sum(s1s1tc1h,2);
                    itidata.s2n.h1h(ii,rpt)=localcoms1h;

                    % TODO one-su com@2h
                    s1tsel2h=ismember(FT_SPIKE.trial{susel},s1trl2h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;

                    s1tfreq2h=histcounts(FT_SPIKE.time{susel}(s1tsel2h)-5-delay,0:0.25:5)./numel(s1trl2h);

                    s2tsel2h=ismember(FT_SPIKE.trial{susel},s2trl2h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;
                    s2tfreq2h=histcounts(FT_SPIKE.time{susel}(s2tsel2h)-5-delay,0:0.25:5)./numel(s2trl2h);

                    s1itimm2h=mean((s1tfreq2h+s2tfreq2h)./2,2);
                    s1s1tc2h=s1tfreq2h-s1itimm2h;
                    s1s1tc2h(s1s1tc2h<0)=0;
                    localcoms2h=sum((1:20).*s1s1tc2h)./sum(s1s1tc2h,2);
                    itidata.s2n.h2h(ii,rpt)=localcoms2h;
                end
            end
            save(fullfile('binary',opt.filename),'itidata');

        end

        function iticorr=corr_stats(itidata)
            % quick dirty stats
            rpt=size(itidata.s1n.h1h,2);
            iticorr=nan(1,rpt);
            for ii=1:rpt
                corrmat=[itidata.s1n.h1h(:,ii),itidata.s1n.h2h(:,ii);...
                    itidata.s2n.h1h(:,ii),itidata.s2n.h2h(:,ii)];
                corr_fini=all(isfinite(corrmat),2);
                r=corr(corrmat(corr_fini,1),corrmat(corr_fini,2));
                iticorr(ii)=r;
            end
        end

        function delay_vs_iti(opt)
            arguments
                opt.filename (1,:) char = 'com_halfs_100.mat'
                opt.plot (1,1) logical = true
                opt.skip_dur (1,1) logical = true
                opt.odor_only (1,1) logical = true
                opt.iti (1,1) logical = false
                opt.iti_filename (1,:) char = 'temp0809_univ_wave.mat'
            end

            load(fullfile('binary',opt.filename),'com_halfs');
            itidata=wave.univ_wave.stats('loadfile',true,'filename',opt.iti_filename)

            
        end


        function nonmemdata=nonmemstats(opt)
            arguments
                opt.rpt (1,1) double =100
                opt.loadfile (1,1) logical = false
                opt.filename (1,:) char = 'temp0809_nonmem.mat'
            end
            error("Not ready")
            if opt.loadfile
                load(fullfile('binary/',opt.filename),'nonmemdata');
                return
            end
            delay=6;
            global_init;
            wrs_mux_meta=ephys.get_wrs_mux_meta();
            % TODO nonmem

            [nonmem_data.s1n.s1iti,nonmem_data.s1n.s2iti,nonmem_data.s2n.s1iti,nonmem_data.s2n.s2iti]=deal([]);

            itidata.s1n.h1h=nan(numel(nonmem_data.s1n.id),100);
            itidata.s1n.h2h=nan(numel(nonmem_data.s1n.id),100);
            itidata.s2n.h1h=nan(numel(nonmem_data.s2n.id),100);
            itidata.s2n.h1h=nan(numel(nonmem_data.s2n.id),100);
            [itidata.s1n.sess,itidata.s1n.cid,itidata.s2n.sess,itidata.s2n.cid]=deal([]);


            sessid=-1;
            for ii=1:numel(nonmem_data.s1n.sess)
                if ~(sessid==nonmem_data.s1n.sess(ii))
                    save(fullfile('binary',opt.filename),'itidata');
                    sessid=nonmem_data.s1n.sess(ii);
                    warning(string(sessid))
                    [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
                    trls=FT_SPIKE.trialinfo;
                    % half-half COM TRLS
                    % ______________________________________________
                    s1trl=find(trls(:,5)==4 & all(trls(:,9:10),2) & trls(:,8)==delay);
                    s2trl=find(trls(:,5)==8 & all(trls(:,9:10),2) & trls(:,8)==delay);
                end

                cid=nonmem_data.s1n.id(ii);
                itidata.s1n.sess=[itidata.s1n.sess;sessid];
                itidata.s1n.cid=[itidata.s1n.cid;cid];
                
                susel=strcmp(FT_SPIKE.label,num2str(cid));
    
                % repeats
                
                for rpt=1:opt.rpt
                    s1trl1h=randsample(s1trl,round(numel(s1trl)./2));
                    s1trl2h=setdiff(s1trl,s1trl1h);

                    s2trl1h=randsample(s2trl,round(numel(s2trl)./2));
                    s2trl2h=setdiff(s2trl,s2trl1h);

                    % one-su com@1h
                    s1tsel1h=ismember(FT_SPIKE.trial{susel},s1trl1h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;

                    s1tfreq1h=histcounts(FT_SPIKE.time{susel}(s1tsel1h)-5-delay,0:0.25:5)./numel(s1trl1h);

                    s2tsel1h=ismember(FT_SPIKE.trial{susel},s2trl1h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;
                    s2tfreq1h=histcounts(FT_SPIKE.time{susel}(s2tsel1h)-5-delay,0:0.25:5)./numel(s2trl1h);

                    s1itimm1h=mean((s1tfreq1h+s2tfreq1h)./2,2);
                    s1s1tc1h=s1tfreq1h-s1itimm1h;
                    s1s1tc1h(s1s1tc1h<0)=0;
                    localcoms1h=sum((1:20).*s1s1tc1h)./sum(s1s1tc1h,2);
                    itidata.s1n.h1h(ii,rpt)=localcoms1h;

                    % TODO one-su com@2h
                    s1tsel2h=ismember(FT_SPIKE.trial{susel},s1trl2h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;

                    s1tfreq2h=histcounts(FT_SPIKE.time{susel}(s1tsel2h)-5-delay,0:0.25:5)./numel(s1trl2h);

                    s2tsel2h=ismember(FT_SPIKE.trial{susel},s2trl2h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;
                    s2tfreq2h=histcounts(FT_SPIKE.time{susel}(s2tsel2h)-5-delay,0:0.25:5)./numel(s2trl2h);

                    s1itimm2h=mean((s1tfreq2h+s2tfreq2h)./2,2);
                    s1s1tc2h=s1tfreq2h-s1itimm2h;
                    s1s1tc2h(s1s1tc2h<0)=0;
                    localcoms2h=sum((1:20).*s1s1tc2h)./sum(s1s1tc2h,2);
                    itidata.s1n.h2h(ii,rpt)=localcoms2h;
                    
                    % if ~all(isfinite([localcoms1h,localcoms2h]),'all')
                    %     keyboard();
                    % end
                end
            end

            sessid=-1;
            for ii=1:numel(nonmem_data.s2n.sess)
                if ~(sessid==nonmem_data.s2n.sess(ii))
                    save(fullfile('binary',opt.filename),'itidata');
                    sessid=nonmem_data.s2n.sess(ii);
                    warning(string(sessid))
                    [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
                    trls=FT_SPIKE.trialinfo;
                    % half-half COM TRLS
                    % ______________________________________________
                    s1trl=find(trls(:,5)==4 & all(trls(:,9:10),2) & trls(:,8)==delay);
                    s2trl=find(trls(:,5)==8 & all(trls(:,9:10),2) & trls(:,8)==delay);
                end

                cid=nonmem_data.s2n.id(ii);
                itidata.s2n.sess=[itidata.s2n.sess;sessid];
                itidata.s2n.cid=[itidata.s2n.cid;cid];

                susel=strcmp(FT_SPIKE.label,num2str(cid));
                % repeats
                
                for rpt=1:opt.rpt
                    s1trl1h=randsample(s1trl,round(numel(s1trl)./2));
                    s1trl2h=setdiff(s1trl,s1trl1h);

                    s2trl1h=randsample(s2trl,round(numel(s2trl)./2));
                    s2trl2h=setdiff(s2trl,s2trl1h);

                    % one-su com@1h
                    s1tsel1h=ismember(FT_SPIKE.trial{susel},s1trl1h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;

                    s1tfreq1h=histcounts(FT_SPIKE.time{susel}(s1tsel1h)-5-delay,0:0.25:5)./numel(s1trl1h);

                    s2tsel1h=ismember(FT_SPIKE.trial{susel},s2trl1h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;
                    s2tfreq1h=histcounts(FT_SPIKE.time{susel}(s2tsel1h)-5-delay,0:0.25:5)./numel(s2trl1h);

                    s1itimm1h=mean((s1tfreq1h+s2tfreq1h)./2,2);
                    s1s1tc1h=s1tfreq1h-s1itimm1h;
                    s1s1tc1h(s1s1tc1h<0)=0;
                    localcoms1h=sum((1:20).*s1s1tc1h)./sum(s1s1tc1h,2);
                    itidata.s2n.h1h(ii,rpt)=localcoms1h;

                    % TODO one-su com@2h
                    s1tsel2h=ismember(FT_SPIKE.trial{susel},s1trl2h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;

                    s1tfreq2h=histcounts(FT_SPIKE.time{susel}(s1tsel2h)-5-delay,0:0.25:5)./numel(s1trl2h);

                    s2tsel2h=ismember(FT_SPIKE.trial{susel},s2trl2h)...
                        & FT_SPIKE.time{susel} >= 5+delay...
                        & FT_SPIKE.time{susel} < 10+delay;
                    s2tfreq2h=histcounts(FT_SPIKE.time{susel}(s2tsel2h)-5-delay,0:0.25:5)./numel(s2trl2h);

                    s1itimm2h=mean((s1tfreq2h+s2tfreq2h)./2,2);
                    s1s1tc2h=s1tfreq2h-s1itimm2h;
                    s1s1tc2h(s1s1tc2h<0)=0;
                    localcoms2h=sum((1:20).*s1s1tc2h)./sum(s1s1tc2h,2);
                    itidata.s2n.h2h(ii,rpt)=localcoms2h;
                end
            end
            save(fullfile('binary',opt.filename),'itidata');

        end
    end
end





% TODO COM curve scale return


