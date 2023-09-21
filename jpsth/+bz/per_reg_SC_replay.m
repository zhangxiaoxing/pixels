load(fullfile('binary','su_meta.mat'));
load(fullfile('binary','wrs_mux_meta.mat'));
load(fullfile('binary','trials_dict.mat'),'trials_dict');
global_init;
[sig,~]=bz.load_sig_sums_conn_file('pair',false);
sig=bz.join_fc_waveid(sig,wrs_mux_meta.wave_id);
congrusel=pct.su_pairs.get_congru(sig.waveid,'odor_only',true);
[gc,gr]=groupcounts(reshape(sig.reg(congrusel,5,:),[],1));
greys=ephys.getGreyRegs('range','grey','mincount',0);
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
greyids=int32(cell2mat(idmap.reg2ccfid.values(greys)));
regsel=all(ismember(sig.reg(:,5,:),intersect(gr(gc>20),greyids)),3);

% scregs=table2cell(cell2table(idmap.ccfid2reg.values(num2cell(intersect(gr(gc>20),greyids)))));

pref_samp=[4 4 8 8 4 8];
np_dur=[6 3 6 3 0 0];

sps=30000;
% per session
for sess=reshape(unique(sig.sess(scids)),1,[])
    disp(sess)
    [spkID,spkTS,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sess,'keep_trial',true,'jagged',true); 
    spkID=uint16(spkID);
    trials=cell2mat(trials_dict(sess));
    session_tick=wave.replay.sessid2length(sess);
    before_sec=(trials(1,1)./sps-60);
    after_sec=(session_tick-trials(end,2))./sps-60-1;
    scids=find(congrusel & regsel & sig.sess==sess);
    for scid=reshape(scids,1,[])
        wids=sig.waveid(scid,:);
        [~,keyidx]=min(wids);
        pref_delay_trls=all(trials(:,9:10)==1,2) & trials(:,5)==pref_samp(wids(keyidx)) & ismember(trials(:,8),setdiff([3 6],np_dur(wids(keyidx))));
        np_delay_trls=all(trials(:,9:10)==1,2) & trials(:,5)==setdiff([4,8],pref_samp(wids(keyidx))) & ismember(trials(:,8),setdiff([3 6],np_dur(wids(keyidx))));

        [prefsec,prefspk,npsec,npspk,pref_iti_sec,pref_iti_spk,np_iti_sec,np_iti_spk]=deal(0);
        for dd=[3 6]
            %pref
            prefsel=pref_delay_trls & trials(:,8)==dd;
            prefsec=prefsec+nnz(prefsel).*dd;
            coupleOne(sig.suid(scid,:),FT_SPIKE,prefsel)
            prefspk=prefspk+nnz(ismember(FT_SPIKE.trial{ftsel},find(prefsel)) & FT_SPIKE.time{ftsel}>=1 & FT_SPIKE.time{ftsel}<4);
            %nonpref
            npsel=np_delay_trls & trials(:,8)==dd;
            npsec=npsec+nnz(npsel).*dd;
            npspk=npspk+nnz(ismember(FT_SPIKE.trial{ftsel},find(npsel)) & FT_SPIKE.time{ftsel}>=1 & FT_SPIKE.time{ftsel}<4);
            %pref iti %
            pref_iti_sec=pref_iti_sec+nnz(prefsel).*9;
            pref_iti_spk=nnz(ismember(FT_SPIKE.trial{ftsel},find(prefsel)) & FT_SPIKE.time{ftsel}>=dd+5 & FT_SPIKE.time{ftsel}<dd+14); % 9s per ITI
            if any(find(prefsel)==size(trials,1)) && (session_tick-trials(end,2))/sps<13
                pref_iti_sec=pref_iti_sec-13+((session_tick-trials(end,2))/sps);
            end
            %np iti
            np_iti_sec=np_iti_sec+nnz(npsel).*9;
            np_iti_spk=nnz(ismember(FT_SPIKE.trial{ftsel},find(npsel)) & FT_SPIKE.time{ftsel}>=dd+5 & FT_SPIKE.time{ftsel}<dd+14); % 9s per ITI
            if any(find(npsel)==size(trials,1)) && (session_tick-trials(end,2))/sps<13
                np_iti_sec=np_iti_sec-13+((session_tick-trials(end,2))/sps);
            end
        end

        %  corresponding network in pre task, post task
        before_sec=trials(1,1)./sps;
        before_spk=nnz(spkID==cid & spkTS<trials(1,1));

        after_sec=(session_tick-trials(end,2))./sps-60;
        if after_sec<=0
            after_sec=nan;
            after_spk=nan;
        else
            after_spk=nnz(spkID==cid & spkTS>=(trials(end,2)+60.*sps));
        end

        cidreg=su_meta.reg_tree{5,su_meta.sess==sess & su_meta.allcid==cid};
        su_sums.(cidreg)=[su_sums.(cidreg);...
            array2table([prefspk,prefsec,npspk,npsec,pref_iti_spk,pref_iti_sec,np_iti_spk,np_iti_sec,before_spk,before_sec,after_spk,after_sec],...
            'VariableNames',{'DelaySPK','DelaySec','NPDelaySPK','NPDelaySec','ITISPK','ITISec','NPITISPK','NPITISec','BeforeSPK','BeforeSec','AfterSPK','AfterSec'})];

    end
end

blame=vcs.blame();
save(fullfile("binary","per_reg_per_su_fr_replay.mat"),"su_sums","blame");

%% plot
fh=figure();
tiledlayout('flow')
for onereg=reshape(statreg,1,[])
    nexttile;
    hold on
    frmat=su_sums.(onereg{1});
    fr=table2array(frmat(:,1:2:end))./table2array(frmat(:,2:2:end));
    mdm=nanmedian(fr);
    ci=bootci(100,@(x) nanmedian(x),fr);
    bar(mdm,'FaceColor','k');
    errorbar(1:6,mdm,ci(1,:)-mdm,ci(2,:)-mdm,'k.');
    set(gca,'XTick',1:6,'XTickLabel',{'Delay','NPDelay','ITI','NPITI','Before','After'});
    ylabel('FR (Hz)');
    title(onereg{1});
end
savefig(fh,fullfile('binary','per_reg_WM_neuron_FR_replay.fig'));


function coupleOne(suid,FT_SPIKE,prefsel,delay)
selL=strcmp(FT_SPIKE.label,num2str(suid(1)));
selF=strcmp(FT_SPIKE.label,num2str(suid(2)));
for t=find(prefsel)
    if delay
        spkL=FT_SPIKE.time{selL}(FT_SPIKE.trial{selL}==t & FT_SPIKE.time{selL}>=1... ...
            & FT_SPIKE.time{selL}<(FT_SPIKE.trialinfo(t,8)+1)
    else %iti

    end

end    
end
