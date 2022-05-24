% TODO rev_fc

%% const
if ~exist('out','var')
    out=struct();
end
sampmap=containers.Map({'s1','s2'},[4,8]);
delay=6;
bin_size=0.25;
[sig,~]=bz.load_sig_pair('CTX',true);
com_map=wave.get_com_map(curve=true,rnd_half=false);
fwd_dir=false;
%% iteration
for sess=1:116
    disp(sess)
    sess_str=sprintf('s%d',sess);
    if ~isfield(com_map,sess_str)
        continue
    end
    for samp=["s1","s2"]
        %% tcom
        sess_suids=sig.suid(sig.sess==sess,:);
        tcom_suid=com_map.(sess_str).(samp).keys();

        lead_follow_sig=[];

        for lead=tcom_suid
            for follow=tcom_suid
                if any(sess_suids(:,1)==lead{1} & sess_suids(:,2)==follow{1})
                    lead_follow_sig=[lead_follow_sig;com_map.(sess_str).(samp)(lead{1}),com_map.(sess_str).(samp)(follow{1}),1];
                else
                    lead_follow_sig=[lead_follow_sig;com_map.(sess_str).(samp)(lead{1}),com_map.(sess_str).(samp)(follow{1}),0];
                end
            end
        end

        %% overall FC figure
        if false
            sigsel=lead_follow_sig(:,3)==1;
            TSspan=max(lead_follow_sig(:,1));
            figure('Color','w')
            hold on
            scatter(lead_follow_sig(~sigsel,1),lead_follow_sig(~sigsel,2),4,[0.7,0.7,0.7],'filled')
            scatter(lead_follow_sig(sigsel,1),lead_follow_sig(sigsel,2),4,'red','filled')
            plot([0,TSspan],[0,TSspan],'--k');
            xlabel('Lead TCOM')
            ylabel('Follow TCOM')
            title(samp)
        end

        %% FCSP freq
        fcsp_suid=intersect(sess_suids(:),cell2mat(tcom_suid));
        [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sess,'keep_trial',true,'suids',unique(fcsp_suid),'only_delay',true);
        %TODO: possible trial number criteria

        fcs=sess_suids(all(ismember(sess_suids,cell2mat(tcom_suid)),2),:);
        fc_count=nan(size(fcs,1),size(trials,1),delay/bin_size);
        for fci=1:size(fcs,1)
            disp(fci)
            for ti=1:size(trials,1)
                susel1=strcmp(FT_SPIKE.label,num2str(fcs(fci,1)));
                susel2=strcmp(FT_SPIKE.label,num2str(fcs(fci,2)));
                for onset_idx=1:size(fc_count,3)
                    onset=(onset_idx-1)*bin_size;
                    spk1=FT_SPIKE.time{susel1}(FT_SPIKE.trial{susel1}==ti & FT_SPIKE.time{susel1}<onset+bin_size & FT_SPIKE.time{susel1}>=onset); %FT_SPIKE.time return [0,6]
                    spk2=FT_SPIKE.time{susel2}(FT_SPIKE.trial{susel2}==ti & FT_SPIKE.time{susel2}<onset+bin_size & FT_SPIKE.time{susel2}>=onset);
                    deltaT=spk2-spk1.'; % of size (n_spk1,n_spk2)
                    if fwd_dir
                        fc_evts=deltaT>=0.0008 & deltaT<0.01; % of size (n_spk1,n_spk2)
                    else
                        fc_evts=deltaT<=-0.0008 & deltaT>-0.01; % of size (n_spk1,n_spk2)
                    end
                    fc_count(fci,ti,onset_idx)=nnz(any(fc_evts)); % counting following spks consistent with fc
                end
            end
        end

        trl_sel=trials(:,5)==sampmap(samp) & trials(:,8)==delay & all(trials(:,9:10),2);
        fc_bin=permute(mean(fc_count(:,trl_sel,:),2),[1,3,2]);
        fc_norm_bin=floor(normalize(fc_bin,2,'range')*255+1);
        fc_tcom=cell2mat(com_map.(sess_str).(samp).values(num2cell(fcs)));
        %% correlation movie
        if false
            figure();
            imagesc(normalize(fc_bin,2,'range'))
            cmap=colormap('turbo');

            fh=figure('Color','w');
            v=VideoWriter('FCSPR_TCOM.mp4','MPEG-4');
            open(v)
            for tt=1:size(fc_norm_bin,2)
                cla
                scatter(fc_tcom(:,1),fc_tcom(:,2),16,fc_norm_bin(:,tt),'filled');
                xline(tt,'--k')
                yline(tt,'--k')
                text(tt,tt,'time','HorizontalAlignment','right','VerticalAlignment','bottom')
                set(gca(),'XTick',0.5:12:24.5,'XTickLabel',0:3:6,'YTick',0.5:12:24.5,'YTickLabel',0:3:6)
                xlim([0.5,24.5])
                ylim([0.5,24.5])
                xlabel('Leading TCOM (s)')
                ylabel('Following TCOM (s)')
                cbh=colorbar();
                cbh.Label.String='Normalized FCSPR';
                set(cbh,'Ticks',0:128:256,'TickLabels',0:0.5:1)
                for i=1:15
                    writeVideo(v,getframe(fh))
                end
            end
            close(v)
        end

        %% curve corr & xcorr
        % lead vs fcspr corr
        % follow vs fcspr corr
        % lead vs follow corr

        corrtype="Pearson";
        corrs=nan(size(fcs,1),3);
        for ii=1:size(fcs,1)
            lead_curve=com_map.(sess_str).(sprintf('%scurve',samp))(fcs(ii,1)).';
            follow_curve=com_map.(sess_str).(sprintf('%scurve',samp))(fcs(ii,2)).';
            fcspr=fc_bin(ii,:).';

            corrs(ii,1)=corr(lead_curve,follow_curve,type=corrtype);
            corrs(ii,2)=corr(lead_curve,fcspr,type=corrtype);
            corrs(ii,3)=corr(follow_curve,fcspr,type=corrtype);
        end
        if fwd_dir
            out.(sess_str).(samp).corrs=corrs;
        else
            out.(sess_str).(samp).rev_corrs=corrs;
        end

        if false
            fh=figure('Color','w');
            bar(mean(corrs))
            set(gca(),'XTick',1:3,'XTickLabel',{'Lead-Follow','Lead-FCSPR','Follow-FCSPR'},'XTickLabelRotation',90)
        end
        % xcorr
        maxlag=4;
        xcorrs=nan(size(fcs,1),3,2*maxlag+1);
        for ii=1:size(fcs,1)
            lead_curve=normalize(com_map.(sess_str).(sprintf('%scurve',samp))(fcs(ii,1)).','range');
            follow_curve=normalize(com_map.(sess_str).(sprintf('%scurve',samp))(fcs(ii,2)).','range');
            fcspr=normalize(fc_bin(ii,:).','range');

            xcorrs(ii,1,:)=xcorr(lead_curve,follow_curve,maxlag);
            xcorrs(ii,2,:)=xcorr(lead_curve,fcspr,maxlag);
            xcorrs(ii,3,:)=xcorr(follow_curve,fcspr,maxlag);
        end
        if fwd_dir
            out.(sess_str).(samp).xcorrs=xcorrs;
        else
            out.(sess_str).(samp).rev_xcorrs=xcorrs;
        end
        if false
            fh=figure('Color','w');
            lbls={'Lead-Follow','Lead-FCSPR','Follow-FCSPR'};
            for subi=1:3
                subplot(1,3,subi)
                hold on
                plot(squeeze(mean(xcorrs(:,subi,:))));
                xline(maxlag+1,'--k')
                title(lbls{subi});
                set(gca,'XTick',1:2*maxlag+1,'XTickLabel',-maxlag:maxlag)
            end
        end
    end
end
% t=rand(1,10);
% a=[t,0,0,0];
% b=[0,0,0,t];
% figure();plot(xcorr(b,a,6))
return
%% all corrs
[xcorrs,corrs,rev_xcorrs,rev_corrs]=deal([]);
sesses=fieldnames(out);
for sess=1:numel(sesses)
    samps=fieldnames(out.(sesses{sess}));
    for samp=1:numel(samps)
        xcorrs=[xcorrs;out.(sesses{sess}).(samps{samp}).xcorrs];
        corrs=[corrs;out.(sesses{sess}).(samps{samp}).corrs];
        
        rev_xcorrs=[rev_xcorrs;out.(sesses{sess}).(samps{samp}).rev_xcorrs];
        rev_corrs=[rev_corrs;out.(sesses{sess}).(samps{samp}).rev_corrs];

    end
end
fh=figure('Color','w');
bar([nanmean(corrs),nanmean(rev_corrs)])
set(gca(),'XTick',1:3,'XTickLabel',{'Lead-Follow','Lead-FCSPR','Follow-FCSPR'},'XTickLabelRotation',90)



fh=figure('Color','w');
lbls={'Lead-Follow','Lead-FCSPR','Follow-FCSPR'};
for subi=1:3
    subplot(1,3,subi)
    hold on
    plot(squeeze(mean(xcorrs(:,subi,:))));
    xline(maxlag+1,'--k')
    title(lbls{subi});
    set(gca,'XTick',1:2*maxlag+1,'XTickLabel',-maxlag:maxlag)
end