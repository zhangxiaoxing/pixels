% TODO: universal wave
delay=6;
global_init;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false,'odor_only',true);
[fh6,imdata]=wave.plot_pct_wave(com_map,'comb_set',4,'flex_sort',true,'scale',[0.1,0.7],'gauss2d',true,'delay',delay);

homedir=ephys.util.getHomedir('type','raw');

[imdata.s1n.s1iti,imdata.s1n.s2iti,imdata.s2n.s1iti,imdata.s2n.s2iti]=deal([]);

sessid=-1;
for ii=1:numel(imdata.s1n.sess)
    if ~(sessid==imdata.s1n.sess(ii))
        sessid=imdata.s1n.sess(ii);
        [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
        trls=FT_SPIKE.trialinfo;
        s1trl=find(trls(:,5)==4 & all(trls(:,9:10),2) & trls(:,8)==delay);
        s2trl=find(trls(:,5)==8 & all(trls(:,9:10),2) & trls(:,8)==delay);
    end
    
    cid=imdata.s1n.id(ii);
    susel=strcmp(FT_SPIKE.label,num2str(cid));
    s1tsel=ismember(FT_SPIKE.trial{susel},s1trl)...
        & FT_SPIKE.time{susel} >= 5+delay...
        & FT_SPIKE.time{susel} < 10+delay;
    s1tfreq=histcounts(FT_SPIKE.time{susel}(s1tsel)-5-delay,0:0.25:5)./numel(s1trl);
    imdata.s1n.s1iti=[imdata.s1n.s1iti;s1tfreq];

    s2tsel=ismember(FT_SPIKE.trial{susel},s2trl)...
        & FT_SPIKE.time{susel} >= 5+delay...
        & FT_SPIKE.time{susel} < 10+delay;
    s2tfreq=histcounts(FT_SPIKE.time{susel}(s2tsel)-5-delay,0:0.25:5)./numel(s2trl);
    imdata.s1n.s2iti=[imdata.s1n.s2iti;s2tfreq];

end

sessid=-1;
for ii=1:numel(imdata.s2n.sess)
    if ~(sessid==imdata.s2n.sess(ii))
        sessid=imdata.s2n.sess(ii);
        [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
        trls=FT_SPIKE.trialinfo;
        s1trl=find(trls(:,5)==4 & all(trls(:,9:10),2) & trls(:,8)==delay);
        s2trl=find(trls(:,5)==8 & all(trls(:,9:10),2) & trls(:,8)==delay);
    end
    
    cid=imdata.s2n.id(ii);
    susel=strcmp(FT_SPIKE.label,num2str(cid));
    s1tsel=ismember(FT_SPIKE.trial{susel},s1trl)...
        & FT_SPIKE.time{susel} >= 5+delay...
        & FT_SPIKE.time{susel} < 10+delay;
    s1tfreq=histcounts(FT_SPIKE.time{susel}(s1tsel)-5-delay,0:0.25:5)./numel(s1trl);
    imdata.s2n.s1iti=[imdata.s2n.s1iti;s1tfreq];

    s2tsel=ismember(FT_SPIKE.trial{susel},s2trl)...
        & FT_SPIKE.time{susel} >= 5+delay...
        & FT_SPIKE.time{susel} < 10+delay;
    s2tfreq=histcounts(FT_SPIKE.time{susel}(s2tsel)-5-delay,0:0.25:5)./numel(s2trl);
    imdata.s2n.s2iti=[imdata.s2n.s2iti;s2tfreq];
end


% TODO COM curve scale return