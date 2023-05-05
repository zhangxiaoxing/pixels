global_init;
wrs_mux_meta=ephys.get_wrs_mux_meta();

[sig,~]=bz.load_sig_sums_conn_file('pair',false);
sig=bz.join_fc_waveid(sig,wrs_mux_meta.wave_id);
sig.congrusel=pct.su_pairs.get_congru(sig.waveid);

usess=unique(sig.sess(sig.congrusel));

sums=[];
for sessid=reshape(usess,1,[])
    disp(sessid)
    usig=sig.suid(sig.sess==sessid & sig.congrusel,:);
    wid=sig.waveid(sig.sess==sessid & sig.congrusel,:);
    [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
    FT_SPIKE.fc_tag=cell(size(FT_SPIKE.timestamp));
    for cidx=1:numel(FT_SPIKE.timestamp)
        FT_SPIKE.fc_tag{cidx}=zeros(size(FT_SPIKE.timestamp{cidx}),'uint8');
    end

    for ss=1:size(usig,1)
        % TODO: pref trials
        leadsel=strcmp(FT_SPIKE.label,num2str(usig(ss,1)));
        folosel=strcmp(FT_SPIKE.label,num2str(usig(ss,2)));
        [pref3,pref6]=bz.rings.preferred_trials_rings(wid(ss,:),trials);
        for tt=reshape([pref3;pref6],1,[])
            leadtssel=FT_SPIKE.trial{leadsel}==tt & FT_SPIKE.time{leadsel}>1 & FT_SPIKE.time{leadsel}<=(trials(tt,8)+1);
            folotssel=FT_SPIKE.trial{folosel}==tt & FT_SPIKE.time{folosel}>1 & FT_SPIKE.time{folosel}<=(trials(tt,8)+1);
            
            % in preferred trial tag
            FT_SPIKE.fc_tag{leadsel}(leadtssel)=bitor(FT_SPIKE.fc_tag{leadsel}(leadtssel),1);
            FT_SPIKE.fc_tag{folosel}(folotssel)=bitor(FT_SPIKE.fc_tag{folosel}(folotssel),1);

            if any(leadtssel) && any(folotssel)
                leadts=FT_SPIKE.timestamp{leadsel}(leadtssel);
                folots=FT_SPIKE.timestamp{folosel}(folotssel);
                latency=folots-leadts.';
                
                leadtag=any(latency>24 & latency<=300,2);
                folotag=any(latency>24 & latency<=300,1);

                if any(leadtag)
                    FT_SPIKE.fc_tag{leadsel}(leadtssel)=bitor(FT_SPIKE.fc_tag{leadsel}(leadtssel),uint8(leadtag*2).');
                end

                if any(folotag)
                    FT_SPIKE.fc_tag{folosel}(folotssel)=bitor(FT_SPIKE.fc_tag{folosel}(folotssel),uint8(folotag*2));
                end
            end
        end
    end

    ufccid=unique(usig);
    for ccid=reshape(ufccid,1,[])
        ccsel=strcmp(FT_SPIKE.label,num2str(ccid));
        onefctag=FT_SPIKE.fc_tag{ccsel};
        sums=[sums;sessid,ccid,nnz(bitand(onefctag,2)),nnz(bitand(onefctag,1)),numel(onefctag)];
    end
end
blame=vcs.blame();
save(fullfile("bzdata","FC_spike_proportion.mat"),'sums','blame')
% prctile(double(sums(:,3))./double(sums(:,4)),[25 50 75])
fc_spk_frac=double(sums(:,3))./double(sums(:,4));
prctile(fc_spk_frac,[25 50 75])

fh=figure();
bh=boxplot(fc_spk_frac,'Colors','k','Whisker',inf,'Widths',0.8);
fh.Children.Children.Children(7).LineStyle='-';
fh.Children.Children.Children(6).LineStyle='-';
% xlim([0.75,1.25])
ylabel('FC spikes / all spikes (%)')
set(gca,'YTick',0:0.2:1,'YTickLabel',0:20:100,'XTick',[])

