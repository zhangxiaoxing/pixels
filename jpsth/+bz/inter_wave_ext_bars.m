function fh=inter_wave_ext_bars(sel_meta)
% [sig,pair]=bz.load_sig_pair('pair',true);
% [sig,pair]=bz.load_sig_sums_conn_file('pair',true);

% meta=ephys.util.load_meta();
% waveid=zeros(size(meta_input.sens_only))+meta_input.sens_only+2.*meta_input.dur_only+4*meta_input.mixed;
% 
% for usess=reshape(unique(meta.sess),1,[])
%     asel=meta.sess==usess;
%     ssel=sig.sess==usess;
%     [~,loc]=ismember(sig.suid(ssel,:),int32(meta.allcid(asel)));
%     sig.waveid(ssel,:)=subsref(waveid(asel),struct(type='()',subs={{loc}}));
% 
%     psel=pair.sess==usess;
%     [~,loc]=ismember(pair.suid(psel,:),int32(meta.allcid(asel)));
%     pair.waveid(psel,:)=subsref(waveid(asel),struct(type='()',subs={{loc}}));
% end

[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
sig=bz.join_fc_waveid(sig,sel_meta.wave_id);
pair=bz.join_fc_waveid(pair,sel_meta.wave_id);


idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
CTX_sel_sig=squeeze(sig.reg(:,2,:)==idmap.reg2ccfid('CTX'));
CNU_sel_sig=squeeze(sig.reg(:,2,:)==idmap.reg2ccfid('CNU'));
BSt_sel_sig=squeeze(sig.reg(:,1,:)==idmap.reg2ccfid('BS'));

CTX_sel_pair=squeeze(pair.reg(:,2,:)==idmap.reg2ccfid('CTX'));
CNU_sel_pair=squeeze(pair.reg(:,2,:)==idmap.reg2ccfid('CNU'));
BSt_sel_pair=squeeze(pair.reg(:,1,:)==idmap.reg2ccfid('BS'));

reg_cell={CTX_sel_sig,CTX_sel_pair;CNU_sel_sig,CNU_sel_pair;BSt_sel_sig,BSt_sel_pair};
%%
[cnt_mix,cnt_dur_only,cnt_sens_only,cnt_nonmem,cnt_ovall,ci_mix,ci_dur_only,ci_sens_only,ci_nonmem,ci_ovall]=deal(nan(3,3,2)); % dimord = from reg, to reg, [sig pair]
[mm_mix,mm_dur_only,mm_sens_only,mm_nonmem,mm_ovall]=deal(nan(3,3));
for from_idx=1:3
    for to_idx=1:3 % [from reg & to reg & waveid=[]],sig then pair
        cnt_mix(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & (all(ismember(sig.waveid,1:4),2))),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & (all(ismember(pair.waveid,1:4),2)))];
        cnt_dur_only(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & (all(ismember(sig.waveid,7:8),2))),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & (all(ismember(pair.waveid,7:8),2)))];
        cnt_sens_only(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & (all(ismember(sig.waveid,5:6),2))),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & (all(ismember(pair.waveid,5:6),2)))];
        cnt_nonmem(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & (all(sig.waveid==0,2))),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & (all(pair.waveid==0,2)))];
        cnt_ovall(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2)),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2))];

        [mm_mix(from_idx,to_idx),ci_mix(from_idx,to_idx,:)]=binofit(cnt_mix(from_idx,to_idx,1),cnt_mix(from_idx,to_idx,2));
        [mm_dur_only(from_idx,to_idx),ci_dur_only(from_idx,to_idx,:)]=binofit(cnt_dur_only(from_idx,to_idx,1),cnt_dur_only(from_idx,to_idx,2));
        [mm_sens_only(from_idx,to_idx),ci_sens_only(from_idx,to_idx,:)]=binofit(cnt_sens_only(from_idx,to_idx,1),cnt_sens_only(from_idx,to_idx,2));
        [mm_nonmem(from_idx,to_idx),ci_nonmem(from_idx,to_idx,:)]=binofit(cnt_nonmem(from_idx,to_idx,1),cnt_nonmem(from_idx,to_idx,2));
        [mm_ovall(from_idx,to_idx),ci_ovall(from_idx,to_idx,:)]=binofit(cnt_ovall(from_idx,to_idx,1),cnt_ovall(from_idx,to_idx,2));
    end
end

%% stats
% sel vs nonmem
ds={cnt_sens_only,cnt_dur_only,cnt_mix};
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
mms={mm_nonmem,mm_sens_only,mm_dur_only,mm_mix};
cis={ci_nonmem,ci_sens_only,ci_dur_only,ci_mix};
fns={'nonmem','sens','dur','mix'};

fh.congru_FC=figure('Color','w','Position',[100,100,1280,500]);
tiledlayout(2,5)
%% stats table
thall=nexttile(5,[2,1]);thall.Visible='off';
th1.Position=thall.Position;th1.Units=thall.Units;th1.Position(4)=thall.Position(4)/3.1;
th2=th1;th2.Position(2)=thall.Position(2)+thall.Position(4)*0.33;
th3=th1;th3.Position(2)=thall.Position(2)+thall.Position(4)*0.67;
ephys.util.figtable(fh.congru_FC,th1,pratio.mix_sel_v_nm,'title','Mixed')
ephys.util.figtable(fh.congru_FC,th2,pratio.dur_sel_v_nm,'title','Duration')
ephys.util.figtable(fh.congru_FC,th3,pratio.sens_sel_v_nm,'title','Sensory')
%%
for ii=1:4
    mm_curr=mms{ii};
    ci_curr=cis{ii};
    nexttile(ii);
    hold on
    bh=bar(mm_curr,'grouped');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,ci_curr(:,1,1).'-bh(1).YEndPoints,ci_curr(:,1,2).'-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,ci_curr(:,2,1).'-bh(2).YEndPoints,ci_curr(:,2,2).'-bh(2).YEndPoints,'k.');
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,ci_curr(:,3,1).'-bh(3).YEndPoints,ci_curr(:,3,2).'-bh(3).YEndPoints,'k.');
    bh(1).FaceColor='w';
    bh(2).FaceColor='r';
    bh(3).FaceColor='b';
    if ii<4
        ylim([0,0.015])
        set(gca(),'YTick',[0,0.005,0.01,0.015],'YTickLabel',[0,0.5,1,1.5],'XTick',[])
    else
        set(gca(),'YTickLabel',get(gca(),'YTick').*100,'XTick',[])
    end
    ylabel('Fucntional coupling rate (%)')
    title(fns{ii});
end

for ii=1:4
    mm_curr=mms{ii}./mm_nonmem;
    ci_curr=cis{ii}./ci_nonmem;
    nexttile(5+ii);
    hold on
    bh=bar(mm_curr,'grouped');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,ci_curr(:,1,1).'-bh(1).YEndPoints,ci_curr(:,1,2).'-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,ci_curr(:,2,1).'-bh(2).YEndPoints,ci_curr(:,2,2).'-bh(2).YEndPoints,'k.');
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,ci_curr(:,3,1).'-bh(3).YEndPoints,ci_curr(:,3,2).'-bh(3).YEndPoints,'k.');
    bh(1).FaceColor='w';
    bh(2).FaceColor='r';
    bh(3).FaceColor='b';
    yline(1,'k--')
    if ii<4
        ylim([0,4])
    end
    set(gca(),'XTick',[])
    ylabel('FC rate relative to non-mem')
    title(fns{ii});
end
sgtitle('grp=from, color=to, LMR|wrb={CTX CNU BST}')
%single modality <-> mixed modality
%%
[cnt_dur2mix,cnt_mix2dur,cnt_sens2mix,cnt_mix2sens,ci_dur2mix,ci_mix2dur,ci_sens2mix,ci_mix2sens]=deal(nan(3,3,2));
[mm_dur2mix,mm_mix2dur,mm_sens2mix,mm_mix2sens]=deal(nan(3,3));
for from_idx=1:3
    for to_idx=1:3
        cnt_dur2mix(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & ismember(sig.waveid(:,1),7:8) & ismember(sig.waveid(:,2),1:4)),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & ismember(pair.waveid(:,1),7:8) & ismember(pair.waveid(:,2),1:4))];
        cnt_mix2dur(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & ismember(sig.waveid(:,1),1:4) & ismember(sig.waveid(:,2),7:8)),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & ismember(pair.waveid(:,1),1:4) & ismember(pair.waveid(:,2),7:8))];
        cnt_sens2mix(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & ismember(sig.waveid(:,1),5:6) & ismember(sig.waveid(:,2),1:4)),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & ismember(pair.waveid(:,1),5:6) & ismember(pair.waveid(:,2),1:4))];
        cnt_mix2sens(from_idx,to_idx,:)=[nnz(reg_cell{from_idx,1}(:,1) & reg_cell{to_idx,1}(:,2) & ismember(sig.waveid(:,1),1:4) & ismember(sig.waveid(:,2),5:6)),...
            nnz(reg_cell{from_idx,2}(:,1) & reg_cell{to_idx,2}(:,2) & ismember(pair.waveid(:,1),1:4) & ismember(pair.waveid(:,2),5:6))];

        [mm_dur2mix(from_idx,to_idx),ci_dur2mix(from_idx,to_idx,:)]=binofit(cnt_dur2mix(from_idx,to_idx,1),cnt_dur2mix(from_idx,to_idx,2));
        [mm_mix2dur(from_idx,to_idx),ci_mix2dur(from_idx,to_idx,:)]=binofit(cnt_mix2dur(from_idx,to_idx,1),cnt_mix2dur(from_idx,to_idx,2));
        [mm_sens2mix(from_idx,to_idx),ci_sens2mix(from_idx,to_idx,:)]=binofit(cnt_sens2mix(from_idx,to_idx,1),cnt_sens2mix(from_idx,to_idx,2));
        [mm_mix2sens(from_idx,to_idx),ci_mix2sens(from_idx,to_idx,:)]=binofit(cnt_mix2sens(from_idx,to_idx,1),cnt_mix2sens(from_idx,to_idx,2));
    end
end

mms={mm_dur2mix,mm_mix2dur,mm_sens2mix,mm_mix2sens};
cis={ci_dur2mix,ci_mix2dur,ci_sens2mix,ci_mix2sens};
fns={'dur2mix','mix2dur','sens2mix','mix2sens'};
fh.single_mix=figure('Color','w','Position',[100,100,1280,500]);
tiledlayout(2,5)
for ii=1:4
    mm_curr=mms{ii};
    ci_curr=cis{ii};
    nexttile(ii);
    hold on
    bh=bar(mm_curr,'grouped');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,ci_curr(:,1,1).'-bh(1).YEndPoints,ci_curr(:,1,2).'-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,ci_curr(:,2,1).'-bh(2).YEndPoints,ci_curr(:,2,2).'-bh(2).YEndPoints,'k.');
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,ci_curr(:,3,1).'-bh(3).YEndPoints,ci_curr(:,3,2).'-bh(3).YEndPoints,'k.');
    bh(1).FaceColor='w';
    bh(2).FaceColor='r';
    bh(3).FaceColor='b';
    ylim([0,0.03])
    set(gca(),'YTick',0:0.01:0.03,'YTickLabel',0:3,'XTick',[])
    ylabel('Fucntional coupling rate (%)')
    title(fns{ii});
end
%% ratio
for ii=1:4
    mm_curr=mms{ii}./mm_nonmem;
    ci_curr=cis{ii}./ci_nonmem;
    nexttile(ii+5);
    hold on
    bh=bar(mm_curr,'grouped');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,ci_curr(:,1,1).'-bh(1).YEndPoints,ci_curr(:,1,2).'-bh(1).YEndPoints,'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,ci_curr(:,2,1).'-bh(2).YEndPoints,ci_curr(:,2,2).'-bh(2).YEndPoints,'k.');
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,ci_curr(:,3,1).'-bh(3).YEndPoints,ci_curr(:,3,2).'-bh(3).YEndPoints,'k.');
    bh(1).FaceColor='w';
    bh(2).FaceColor='r';
    bh(3).FaceColor='b';
    yline(1,'k--')
    ylim([0,6])
    set(gca(),'XTick',[])
    ylabel('FC rate relative to non-mem')
    title(fns{ii});
end

%% stats
% sel vs nonmem
ds={cnt_dur2mix,cnt_mix2dur,cnt_sens2mix,cnt_mix2sens};
[pratio.dur2mix,pratio.mix2dur,pratio.sens2mix,pratio.mix2sens]=deal(nan(3,3));
ps={'dur2mix','mix2dur','sens2mix','mix2sens'};

for cnt_idx=1:4
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

%% stats table
thall=nexttile(5,[2,1]);thall.Visible='off';
th1.Position=thall.Position;th1.Units=thall.Units;th1.Position(4)=thall.Position(4)/4.1;
th2=th1;th2.Position(2)=thall.Position(2)+thall.Position(4)*0.25;
th3=th1;th3.Position(2)=thall.Position(2)+thall.Position(4)*0.5;
th4=th1;th4.Position(2)=thall.Position(2)+thall.Position(4)*0.75;
ephys.util.figtable(fh.single_mix,th1,pratio.mix2dur,'title','Mix2Dur')
ephys.util.figtable(fh.single_mix,th2,pratio.dur2mix,'title','Dur2Mix')
ephys.util.figtable(fh.single_mix,th3,pratio.mix2sens,'title','Mix2Sens')
ephys.util.figtable(fh.single_mix,th4,pratio.sens2mix,'title','Sens2Mix')

end

