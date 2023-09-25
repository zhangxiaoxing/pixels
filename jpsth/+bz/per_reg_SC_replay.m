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
regs=cellfun(@(x) idmap.ccfid2reg(x),num2cell(intersect(gr(gc>20),greyids)));
scregs.from=cell2struct(cell(numel(regs),1),regs,1);
scregs.to=cell2struct(cell(numel(regs),1),regs,1);

SPratio.lead=cell2struct(cell(numel(regs),1),regs,1);
SPratio.follow=cell2struct(cell(numel(regs),1),regs,1);

pref_samp=[4 4 8 8 4 8];
np_dur=[6 3 6 3 0 0];

sps=30000;
% per session
reg_pairs=cell(0,2);
for sess=reshape(unique(sig.sess(congrusel & regsel)),1,[])
    disp(sess)

    [spkID,spkTS,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sess,'keep_trial',true,'jagged',true); 
    spkID=uint16(spkID);

    trials=cell2mat(trials_dict(sess));
    session_tick=wave.replay.sessid2length(sess);
    before_sec=(trials(1,1)./sps-60);
    after_sec=(session_tick-trials(end,2))./sps-60-1;
    scids=find(congrusel & regsel & sig.sess==sess);
    for scid=reshape(scids,1,[])
        screg=cellfun(@(x)idmap.ccfid2reg(x),num2cell(squeeze(sig.reg(scid,5,:))));
        reg_pairs=[reg_pairs;screg.'];
        wids=sig.waveid(scid,:);
        [~,keyidx]=min(wids);
        pref_delay_trls=all(trials(:,9:10)==1,2) & trials(:,5)==pref_samp(wids(keyidx)) & ismember(trials(:,8),setdiff([3 6],np_dur(wids(keyidx))));
        np_delay_trls=all(trials(:,9:10)==1,2) & trials(:,5)==setdiff([4,8],pref_samp(wids(keyidx))) & ismember(trials(:,8),setdiff([3 6],np_dur(wids(keyidx))));

        %pref
        prefsec=sum(trials(pref_delay_trls,8));
        [prefSP,pref_lead,pref_follow]=coupleOne(sig.suid(scid,:),FT_SPIKE,pref_delay_trls,true);
        %nonpref
        npsec=sum(trials(np_delay_trls,8));
        [npSP,np_lead,np_follow]=coupleOne(sig.suid(scid,:),FT_SPIKE,np_delay_trls,true);
        %pref iti %
        pref_iti_sec=nnz(pref_delay_trls).*9;
        [pref_iti_SP,iti_lead,iti_follow]=coupleOne(sig.suid(scid,:),FT_SPIKE,pref_delay_trls,false);

        if any(find(pref_delay_trls)==size(trials,1)) && (session_tick-trials(end,2))/sps<13
            pref_iti_sec=pref_iti_sec-13+((session_tick-trials(end,2))/sps);
        end
        %np iti
        np_iti_sec=nnz(np_delay_trls).*9;
        [np_iti_SP,npiti_lead,npiti_follow]=coupleOne(sig.suid(scid,:),FT_SPIKE,np_delay_trls,false);

        if any(find(np_delay_trls)==size(trials,1)) && (session_tick-trials(end,2))/sps<13
            np_iti_sec=np_iti_sec-13+((session_tick-trials(end,2))/sps);
        end

        %  corresponding network in pre task, post task
        before_sec=trials(1,1)./sps;
        [before_SP,before_lead,before_follow]=coupleOuttask(sig.suid(scid,:),spkID,spkTS,0,trials(1,1));

        after_sec=(session_tick-trials(end,2))./sps-60;
        if after_sec<=0
            after_sec=nan;
            after_SP=nan;
        else
            [after_SP,after_lead,after_follow]=coupleOuttask(sig.suid(scid,:),spkID,spkTS,trials(end,2)+60.*sps,session_tick);
        end
        %TODO gather result here

        

        scregs.from.(screg{1})=[scregs.from.(screg{1});...
                        array2table([prefSP,prefsec,npSP,npsec,pref_iti_SP,pref_iti_sec,np_iti_SP,np_iti_sec,before_SP,before_sec,after_SP,after_sec],...
            'VariableNames',{'DelaySP','DelaySec','NPDelaySP','NPDelaySec','ITISP','ITISec','NPITISP','NPITISec','BeforeSP','BeforeSec','AfterSP','AfterSec'})];

        scregs.to.(screg{2})=[scregs.to.(screg{2});...
                        array2table([prefSP,prefsec,npSP,npsec,pref_iti_SP,pref_iti_sec,np_iti_SP,np_iti_sec,before_SP,before_sec,after_SP,after_sec],...
            'VariableNames',{'DelaySP','DelaySec','NPDelaySP','NPDelaySec','ITISP','ITISec','NPITISP','NPITISec','BeforeSP','BeforeSec','AfterSP','AfterSec'})];

        SPratio.lead.(screg{1})=[SPratio.lead.(screg{1});...
                                    array2table([prefSP,pref_lead,npSP,np_lead,pref_iti_SP,iti_lead,np_iti_SP,npiti_lead,before_SP,before_lead,after_SP,after_lead],...
            'VariableNames',{'DelaySP','DelaySPK','NPDelaySP','NPDelaySPK','ITISP','ITISPK','NPITISP','NPITISPK','BeforeSP','BeforeSPK','AfterSP','AfterSPK'})];

        SPratio.follow.(screg{2})=[SPratio.follow.(screg{2});...
            array2table([prefSP,pref_follow,npSP,np_follow,pref_iti_SP,iti_follow,np_iti_SP,npiti_follow,before_SP,before_follow,after_SP,after_follow],...
            'VariableNames',{'DelaySP','DelaySPK','NPDelaySP','NPDelaySPK','ITISP','ITISPK','NPITISP','NPITISPK','BeforeSP','BeforeSPK','AfterSP','AfterSPK'})];

    end
end

blame=vcs.blame();
save(fullfile("binary","per_reg_SC_replay.mat"),"scregs","SPratio","reg_pairs","blame");

%% plot
for dir=["from","to"]
    fh.(dir)=figure();
    tiledlayout('flow')
    for onereg=reshape(sort(regs),1,[])
        nexttile;
        hold on
        frmat=scregs.(dir).(onereg{1});
        fr=table2array(frmat(:,1:2:end))./table2array(frmat(:,2:2:end));
        mdm=nanmedian(fr);
        ci=bootci(100,@(x) nanmedian(x),fr);
        bar(mdm,'FaceColor','k');
        errorbar(1:6,mdm,ci(1,:)-mdm,ci(2,:)-mdm,'k.');
        set(gca,'XTick',1:6,'XTickLabel',{'Delay','NPDelay','ITI','NPITI','Before','After'});
        ylabel('FR (Hz)');
        title(onereg{1});
    end
    sgtitle("SCs "+dir)
end
savefig(fh.from,fullfile('binary','per_reg_SC_from_replay.fig'));
savefig(fh.to,fullfile('binary','per_reg_SC_to_replay.fig'));



%% plot
for dir=["lead","follow"]
    fh.(dir)=figure();
    tiledlayout('flow')
    for onereg=reshape(sort(regs),1,[])
        nexttile;
        hold on
        frmat=SPratio.(dir).(onereg{1});
        fr=table2array(frmat(:,1:2:end))./table2array(frmat(:,2:2:end));
        mdm=nanmedian(fr);
        ci=bootci(100,@(x) nanmedian(x),fr);
        bar(mdm,'FaceColor','k');
        errorbar(1:6,mdm,ci(1,:)-mdm,ci(2,:)-mdm,'k.');
        set(gca,'XTick',1:6,'XTickLabel',{'Delay','NPDelay','ITI','NPITI','Before','After'});
        ylabel('Ratio');
        title(onereg{1});
    end
    sgtitle("SP/all "+dir+" spikes ")
end
savefig(fh.lead,fullfile('binary','per_reg_SC_ratio_lead_replay.fig'));
savefig(fh.follow,fullfile('binary','per_reg_SC_ratio_follow_replay.fig'));





function [SPcnt,leadcnt,followcnt]=coupleOne(suid,FT_SPIKE,prefsel,delay)
selL=strcmp(FT_SPIKE.label,num2str(suid(1)));
selF=strcmp(FT_SPIKE.label,num2str(suid(2)));
[SPcnt,leadcnt,followcnt]=deal(0);
for t=reshape(find(prefsel),1,[])
    if delay
        spkL=FT_SPIKE.time{selL}(FT_SPIKE.trial{selL}==t & FT_SPIKE.time{selL}>=1 ...
            & FT_SPIKE.time{selL}<(FT_SPIKE.trialinfo(t,8)+1));
        spkF=FT_SPIKE.time{selF}(FT_SPIKE.trial{selF}==t & FT_SPIKE.time{selF}>=1 ...
            & FT_SPIKE.time{selF}<(FT_SPIKE.trialinfo(t,8)+1));
    else %iti
        spkL=FT_SPIKE.time{selL}(FT_SPIKE.trial{selL}==t & FT_SPIKE.time{selL}>=(FT_SPIKE.trialinfo(t,8)+5) ...
            & FT_SPIKE.time{selL}<(FT_SPIKE.trialinfo(t,8)+14));
        spkF=FT_SPIKE.time{selF}(FT_SPIKE.trial{selF}==t & FT_SPIKE.time{selF}>=(FT_SPIKE.trialinfo(t,8)+5) ...
            & FT_SPIKE.time{selF}<(FT_SPIKE.trialinfo(t,8)+14));
    end
    lagmat=(spkL-spkF.');
    SPcnt=SPcnt+nnz(lagmat>0 & lagmat<0.01);
    leadcnt=leadcnt+numel(spkL);
    followcnt=followcnt+numel(spkF);
end
end


function [SPcnt,leadcnt,followcnt]=coupleOuttask(suid,spkID,spkTS,onset,offset)
SPcnt=0;
Lts=spkTS(spkID==suid(1) & spkTS>=onset & spkTS<offset);
Fts=spkTS(spkID==suid(2) & spkTS>=onset & spkTS<offset);

leadcnt=numel(Lts);
followcnt=numel(Fts);

if offset<=600000
    Lbin=ones(size(Lts));
    Fbin=ones(size(Fts));
else
    [~,~,Lbin]=histcounts(Lts,0:600000:offset); %20s bin
    [~,~,Fbin]=histcounts(Fts,0:600000:offset); %20s bin
end
for b=reshape(intersect(Lbin,Fbin),1,[])
    ltsb=Lts(Lbin==b);
    ftsb=Fts(Fbin==b);
    lagmat=(ltsb-ftsb.');
    SPcnt=SPcnt+nnz(lagmat>0 & lagmat<300);
end
end

