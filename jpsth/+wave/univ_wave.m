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


s1itimm=mean((imdata.s1n.s1iti(:,5:28)+imdata.s1n.s2iti(:,5:28))./2,2);
S=max(imdata.s1n.s1iti(:,5:28)-s1itimm,[],2);
s1s1t=(imdata.s1n.s1iti-s1itimm)./S;
s1s2t=(imdata.s1n.s2iti-s1itimm)./S;

s2itimm=mean((imdata.s2n.s1iti(:,5:28)+imdata.s2n.s2iti(:,5:28))./2,2);
S=max(imdata.s2n.s1iti(:,5:28)-s2itimm,[],2);
s2s1t=(imdata.s2n.s1iti-s2itimm)./S;
s2s2t=(imdata.s2n.s2iti-s2itimm)./S;



immat={s1s1t(imdata.s1n.sorts,:),s1s2t(imdata.s1n.sorts,:),...
     s2s1t(imdata.s2n.sorts,:),s2s2t(imdata.s2n.sorts,:)};

%% PLOT
scale=[0,0.8];
gk = fspecial('gaussian', [9 3], 1);
figure()
tiledlayout(2,2)
for datum=immat
    nexttile()
    imagesc(conv2(datum{1},gk,'same'),scale)
    set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
    xlim([0.5,size(datum{1},2)+0.5])
    ylim([0.5,size(datum{1},1)+0.5])
    set(gca(),'YDir','normal')
end
colormap("turbo")





% TODO COM curve scale return