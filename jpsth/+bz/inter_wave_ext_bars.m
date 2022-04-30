function fh=inter_wave_ext_bars(anovameta)

[sig,pair]=bz.load_sig_pair('pair',true);

meta=ephys.util.load_meta();
% [dur_sense_mix,dur_exclu,~,sens_exclu]=ephys.get_dul_sel();
% waveid=zeros(size(dur_sense_mix))+4.*dur_sense_mix+2.*dur_exclu+sens_exclu;
waveid=zeros(size(anovameta.sens_only))+4.*anovameta.sens_only+2.*anovameta.dur_only+anovameta.mixed;

for usess=reshape(unique(meta.sess),1,[])
    asel=meta.sess==usess;
    ssel=sig.sess==usess;
    %TODO: update wave component id with dur_sel routine
    [~,loc]=ismember(sig.suid(ssel,:),int32(meta.allcid(asel)));
    sig.wave_id(ssel,:)=subsref(waveid(asel),struct(type='()',subs={{loc}}));

    psel=pair.sess==usess;
    [~,loc]=ismember(pair.suid(psel,:),int32(meta.allcid(asel)));
    pair.wave_id(psel,:)=subsref(waveid(asel),struct(type='()',subs={{loc}}));
end

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
CTX_sel_sig=squeeze(sig.reg(:,2,:)==idmap.reg2ccfid('CTX'));
CNU_sel_sig=squeeze(sig.reg(:,2,:)==idmap.reg2ccfid('CNU'));
BSt_sel_sig=squeeze(sig.reg(:,1,:)==idmap.reg2ccfid('BS'));

CTX_sel_pair=squeeze(pair.reg(:,2,:)==idmap.reg2ccfid('CTX'));
CNU_sel_pair=squeeze(pair.reg(:,2,:)==idmap.reg2ccfid('CNU'));
BSt_sel_pair=squeeze(pair.reg(:,1,:)==idmap.reg2ccfid('BS'));


reg_cell={CTX_sel_sig,CTX_sel_pair;CNU_sel_sig,CNU_sel_pair;BSt_sel_sig,BSt_sel_pair};
% ttl={'Different region','Same region','Overall'};
%%
[cnt_both,cnt_dur_only,cnt_sens_only,cnt_nonmem,cnt_ovall,ci_both,ci_dur_only,ci_sens_only,ci_nonmem,ci_ovall]=deal(nan(3,3,2));
[mm_both,mm_dur_only,mm_sens_only,mm_nonmem,mm_ovall]=deal(nan(3,3));
for from_idx=1:3
    for to_idx=1:3
        cnt_both(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & (all(sig.wave_id==4,2))),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & (all(pair.wave_id==4,2)))];
        cnt_dur_only(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & (all(sig.wave_id==2,2))),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & (all(pair.wave_id==2,2)))];
        cnt_sens_only(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & (all(sig.wave_id==1,2))),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & (all(pair.wave_id==1,2)))];
        cnt_nonmem(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & (all(sig.wave_id==0,2))),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & (all(pair.wave_id==0,2)))];
        cnt_ovall(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2)),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2))];

        [mm_both(from_idx,to_idx),ci_both(from_idx,to_idx,:)]=binofit(cnt_both(from_idx,to_idx,1),cnt_both(from_idx,to_idx,2));
        [mm_dur_only(from_idx,to_idx),ci_dur_only(from_idx,to_idx,:)]=binofit(cnt_dur_only(from_idx,to_idx,1),cnt_dur_only(from_idx,to_idx,2));
        [mm_sens_only(from_idx,to_idx),ci_sens_only(from_idx,to_idx,:)]=binofit(cnt_sens_only(from_idx,to_idx,1),cnt_sens_only(from_idx,to_idx,2));
        [mm_nonmem(from_idx,to_idx),ci_nonmem(from_idx,to_idx,:)]=binofit(cnt_nonmem(from_idx,to_idx,1),cnt_nonmem(from_idx,to_idx,2));
        [mm_ovall(from_idx,to_idx),ci_ovall(from_idx,to_idx,:)]=binofit(cnt_ovall(from_idx,to_idx,1),cnt_ovall(from_idx,to_idx,2));
    end
end

%% stats
% sel vs nonmem
ds={cnt_sens_only,cnt_dur_only,cnt_both};
[pratio.sens_sel_v_nm,pratio.dur_sel_v_nm,pratio.mix_sel_v_nm]=deal(nan(3,3));
ps={'sens_sel_v_nm','dur_sel_v_nm','mix_sel_v_nm'};

for cnt_idx=1:3
    for from_idx=1:3
        for to_idx=1:3
            nm_d=cnt_nonmem(from_idx,to_idx,:);
            sel_d=ds{cnt_idx}(from_idx,to_idx,:);
            [~,~,p]=crosstab([zeros(1,nm_d(2)),ones(1,sel_d(2))],[1:nm_d(2)>nm_d(1),1:sel_d(2)>sel_d(1)]);
%             disp([from_idx,to_idx,p]);
            pratio.(ps{cnt_idx})(from_idx,to_idx)=p;
        end
    end
end

%%


mms={mm_nonmem,mm_sens_only,mm_dur_only,mm_both};
cis={ci_nonmem,ci_sens_only,ci_dur_only,ci_both};
fns={'nonmem','sens','dur','both'};
fh.congru_FC_rate_raw=figure('Color','w','Position',[100,100,1000,205]);
for ii=1:4
    mm_curr=mms{ii};
    ci_curr=cis{ii};
    subplot(1,4,ii);
    hold on
    bh=bar(mm_curr,'grouped');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,ci_curr(:,1,1).'-bh(1).YEndPoints,ci_curr(:,1,2).'-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,ci_curr(:,2,1).'-bh(2).YEndPoints,ci_curr(:,2,2).'-bh(2).YEndPoints,'k.');
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,ci_curr(:,3,1).'-bh(3).YEndPoints,ci_curr(:,3,2).'-bh(3).YEndPoints,'k.');
    bh(1).FaceColor='w';
    bh(2).FaceColor='r';
    bh(3).FaceColor='b';
%     ylim([0,0.015])
    set(gca(),'YTick',[0,0.005,0.01,0.015],'YTickLabel',[0,0.5,1,1.5],'XTick',[])
    ylabel('Fucntional coupling rate (%)')
    title(fns{ii});
end
exportgraphics(fh.congru_FC_rate_raw,'FC_CTX_CNU_BS_raw.pdf','ContentType','vector','Append',true)
fh.congru_FC_ratio=figure('Color','w','Position',[100,100,1000,205]);
for ii=1:4
    mm_curr=mms{ii}./mm_nonmem;
    ci_curr=cis{ii}./ci_nonmem;
    subplot(1,4,ii);
    hold on
    bh=bar(mm_curr,'grouped');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,ci_curr(:,1,1).'-bh(1).YEndPoints,ci_curr(:,1,2).'-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,ci_curr(:,2,1).'-bh(2).YEndPoints,ci_curr(:,2,2).'-bh(2).YEndPoints,'k.');
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,ci_curr(:,3,1).'-bh(3).YEndPoints,ci_curr(:,3,2).'-bh(3).YEndPoints,'k.');
    bh(1).FaceColor='w';
    bh(2).FaceColor='r';
    bh(3).FaceColor='b';
    yline(1,'k--')
%     ylim([0,3])
    set(gca(),'XTick',[])
    ylabel('FC rate relative to non-mem')
    title(fns{ii});
end
exportgraphics(fh.congru_FC_ratio,'FC_CTX_CNU_BS_raw.pdf','ContentType','vector','Append',true)
%single modality <-> mixed modality
%%
[cnt_dur2both,cnt_both2dur,cnt_sens2both,cnt_both2sens,ci_dur2both,ci_both2dur,ci_sens2both,ci_both2sens]=deal(nan(3,3,2));
[mm_dur2both,mm_both2dur,mm_sens2both,mm_both2sens]=deal(nan(3,3));
for from_idx=1:3
    for to_idx=1:3
        cnt_dur2both(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & sig.wave_id(:,1)==2 & sig.wave_id(:,2)==4),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & pair.wave_id(:,1)==2 & pair.wave_id(:,2)==4)];
        cnt_both2dur(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & sig.wave_id(:,1)==4 & sig.wave_id(:,2)==2),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & pair.wave_id(:,1)==4 & pair.wave_id(:,2)==2)];
        cnt_sens2both(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & sig.wave_id(:,1)==1 & sig.wave_id(:,2)==4),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & pair.wave_id(:,1)==1 & pair.wave_id(:,2)==4)];
        cnt_both2sens(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & sig.wave_id(:,1)==4 & sig.wave_id(:,2)==1),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & pair.wave_id(:,1)==4 & pair.wave_id(:,2)==1)];

        [mm_dur2both(from_idx,to_idx),ci_dur2both(from_idx,to_idx,:)]=binofit(cnt_dur2both(from_idx,to_idx,1),cnt_dur2both(from_idx,to_idx,2));
        [mm_both2dur(from_idx,to_idx),ci_both2dur(from_idx,to_idx,:)]=binofit(cnt_both2dur(from_idx,to_idx,1),cnt_both2dur(from_idx,to_idx,2));
        [mm_sens2both(from_idx,to_idx),ci_sens2both(from_idx,to_idx,:)]=binofit(cnt_sens2both(from_idx,to_idx,1),cnt_sens2both(from_idx,to_idx,2));
        [mm_both2sens(from_idx,to_idx),ci_both2sens(from_idx,to_idx,:)]=binofit(cnt_both2sens(from_idx,to_idx,1),cnt_both2sens(from_idx,to_idx,2));
    end
end

mms={mm_dur2both,mm_both2dur,mm_sens2both,mm_both2sens};
cis={ci_dur2both,ci_both2dur,ci_sens2both,ci_both2sens};
fns={'dur2both','both2dur','sens2both','both2sens'};
fh.single_mix_raw=figure('Color','w','Position',[100,100,1000,205]);
for ii=1:4
    mm_curr=mms{ii};
    ci_curr=cis{ii};
    subplot(1,4,ii);
    hold on
    bh=bar(mm_curr,'grouped');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,ci_curr(:,1,1).'-bh(1).YEndPoints,ci_curr(:,1,2).'-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,ci_curr(:,2,1).'-bh(2).YEndPoints,ci_curr(:,2,2).'-bh(2).YEndPoints,'k.');
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,ci_curr(:,3,1).'-bh(3).YEndPoints,ci_curr(:,3,2).'-bh(3).YEndPoints,'k.');
    bh(1).FaceColor='w';
    bh(2).FaceColor='r';
    bh(3).FaceColor='b';
%     ylim([0,0.015])
    set(gca(),'YTick',[0,0.005,0.01,0.015],'YTickLabel',[0,0.5,1,1.5],'XTick',[])
    ylabel('Fucntional coupling rate (%)')
    title(fns{ii});
end
exportgraphics(fh.single_mix_raw,'FC_CTX_CNU_BS_raw.pdf','ContentType','vector','Append',true)

fh.single_mix_ratio=figure('Color','w','Position',[100,100,1000,205]);
for ii=1:4
    mm_curr=mms{ii}./mm_nonmem;
    ci_curr=cis{ii}./ci_nonmem;
    subplot(1,4,ii);
    hold on
    bh=bar(mm_curr,'grouped');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,ci_curr(:,1,1).'-bh(1).YEndPoints,ci_curr(:,1,2).'-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,ci_curr(:,2,1).'-bh(2).YEndPoints,ci_curr(:,2,2).'-bh(2).YEndPoints,'k.');
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,ci_curr(:,3,1).'-bh(3).YEndPoints,ci_curr(:,3,2).'-bh(3).YEndPoints,'k.');
    bh(1).FaceColor='w';
    bh(2).FaceColor='r';
    bh(3).FaceColor='b';
    yline(1,'k--')
%     ylim([0,3])
    set(gca(),'XTick',[])
    ylabel('FC rate relative to non-mem')
    title(fns{ii});
end
exportgraphics(fh.single_mix_ratio,'FC_CTX_CNU_BS_raw.pdf','ContentType','vector','Append',true)
end

