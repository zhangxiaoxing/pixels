early_bins=true;
%% 3s preferred, 6s preferred
dursel=wave.delay_dur_selective('window',1:2);
%[1:fid,2:sess,3:suid,4-6:s1p,s2p,s1|s2p,7-8:FR3,FR6,9:sel_idx,10:AUC
dursel(:,11)=0;
dursel(dursel(:,6)<0.05 & dursel(:,7)>dursel(:,8),11)=3;
dursel(dursel(:,6)<0.05 & dursel(:,7)<dursel(:,8),11)=6;
% to map
dur_sel_map=containers.Map(num2cell(uint64(dursel(:,2)*100000)+uint64(dursel(:,3))),num2cell(dursel(:,11)));


%% 3s selective, 6s selective

meta3=ephys.util.load_meta('delay',3);
meta6=ephys.util.load_meta('delay',6);
if early_bins
    wave_3_6=zeros(size(meta3.sess));
    wave_3_6(((ismember(meta3.mem_type,1:2) & ismember(meta6.mem_type,1:2))...
        |(ismember(meta3.mem_type,3:4) & ismember(meta6.mem_type,3:4)))...
        & any(meta6.per_bin(1:2,:)) & any(meta3.per_bin(1:2,:)))=3;
    wave_3_6(meta3.mem_type>0 & meta6.mem_type==0 & any(meta3.per_bin(1:2,:)))=1;
    wave_3_6(meta3.mem_type==0 & meta6.mem_type>0 & any(meta6.per_bin(1:2,:)))=2;
else
    wave_3_6=zeros(size(meta3.sess));
    wave_3_6(((ismember(meta3.mem_type,1:2) & ismember(meta6.mem_type,1:2))...
        |(ismember(meta3.mem_type,3:4) & ismember(meta6.mem_type,3:4)))...
        )=3;
    wave_3_6(meta3.mem_type>0 & meta6.mem_type==0)=1;
    wave_3_6(meta3.mem_type==0 & meta6.mem_type>0)=2;
end
% to map
wave_id_map=containers.Map(num2cell(...
    uint64(meta3.sess*100000)+uint64(meta3.allcid)),...
    num2cell(wave_3_6));

%% functional coupling
[sig,pair]=bz.load_sig_pair('pair',true);

%sig 3-3 pair 3-3; 3-6; 6-3; 6-6

sig_dur=[cell2mat(dur_sel_map.values(num2cell(int64(sig.sess)*100000+int64(sig.suid(:,1))))),...
    cell2mat(dur_sel_map.values(num2cell(int64(sig.sess)*100000+int64(sig.suid(:,2)))))];

pair_dur=[cell2mat(dur_sel_map.values(num2cell(int64(pair.sess)*100000+int64(pair.suid(:,1))))),...
    cell2mat(dur_sel_map.values(num2cell(int64(pair.sess)*100000+int64(pair.suid(:,2)))))];

sig_wave=[cell2mat(wave_id_map.values(num2cell(int64(sig.sess)*100000+int64(sig.suid(:,1))))),...
    cell2mat(wave_id_map.values(num2cell(int64(sig.sess)*100000+int64(sig.suid(:,2)))))];

pair_wave=[cell2mat(wave_id_map.values(num2cell(int64(pair.sess)*100000+int64(pair.suid(:,1))))),...
    cell2mat(wave_id_map.values(num2cell(int64(pair.sess)*100000+int64(pair.suid(:,2)))))];

% [phat33,pci33]=binofit(nnz(sig_dur(:,1)==3 & sig_wave(:,2)==1)+nnz(sig_dur(:,2)==3 & sig_wave(:,1)==1),...
%     nnz(pair_dur(:,1)==3 & pair_wave(:,2)==1)+nnz(pair_dur(:,2)==3 & pair_wave(:,1)==1));

% [phat3b,pci3b]=binofit(nnz(sig_dur(:,1)==3 & sig_wave(:,2)==3)+nnz(sig_dur(:,2)==3 & sig_wave(:,1)==3),...
%     nnz(pair_dur(:,1)==3 & pair_wave(:,2)==3)+nnz(pair_dur(:,2)==3 & pair_wave(:,1)==3));

% [phat66,pci66]=binofit(nnz(sig_dur(:,1)==6 & sig_wave(:,2)==2)+nnz(sig_dur(:,2)==6 & sig_wave(:,1)==2),...
%     nnz(pair_dur(:,1)==6 & pair_wave(:,2)==2)+nnz(pair_dur(:,2)==6 & pair_wave(:,1)==2));

% [phat6b,pci6b]=binofit(nnz(sig_dur(:,1)==6 & sig_wave(:,2)==3)+nnz(sig_dur(:,2)==6 & sig_wave(:,1)==3),...
%     nnz(pair_dur(:,1)==6 & pair_wave(:,2)==3)+nnz(pair_dur(:,2)==6 & pair_wave(:,1)==3));


nnz_sig_33b=nnz(sig_dur(:,1)==3 & (ismember(sig_wave(:,2),[1 3])))+nnz(sig_dur(:,2)==3 & (ismember(sig_wave(:,1),[1 3])));
nnz_pair_33b=nnz(pair_dur(:,1)==3 & (ismember(pair_wave(:,2),[1 3])))+nnz(pair_dur(:,2)==3 & (ismember(pair_wave(:,1),[1 3])));
nnz_sig_36=nnz(sig_dur(:,1)==3 & sig_wave(:,2)==2)+nnz(sig_dur(:,2)==3 & sig_wave(:,1)==2);
nnz_pair_36=nnz(pair_dur(:,1)==3 & pair_wave(:,2)==2)+nnz(pair_dur(:,2)==3 & pair_wave(:,1)==2);

nnz_sig_63=nnz(sig_dur(:,1)==6 & sig_wave(:,2)==1)+nnz(sig_dur(:,2)==6 & sig_wave(:,1)==1);
nnz_sig_66b=nnz(sig_dur(:,1)==6 & (ismember(sig_wave(:,2),[2 3])))+nnz(sig_dur(:,2)==6 & (ismember(sig_wave(:,1),[2 3])));
nnz_pair_63=nnz(pair_dur(:,1)==6 & pair_wave(:,2)==1)+nnz(pair_dur(:,2)==6 & pair_wave(:,1)==1);
nnz_pair_66b=nnz(pair_dur(:,1)==6 & (ismember(pair_wave(:,2),[2 3])))+nnz(pair_dur(:,2)==6 & (ismember(pair_wave(:,1),[2 3])));

[phat33b,pci33b]=binofit(nnz_sig_33b,nnz_pair_33b);
[phat36,pci36]=binofit(nnz_sig_36,nnz_pair_36);

[phat63,pci63]=binofit(nnz_sig_63,nnz_pair_63);
[phat66b,pci66b]=binofit(nnz_sig_66b,nnz_pair_66b);


fh=figure('Color','w','Position',[32,32,392,417]);
hold on
mm=[phat33b,phat36,phat63,phat66b].*100;
bh=bar([1,2,4,5],diag(mm),'stacked','FaceColor','w');
[bh(2).FaceColor,bh(3).FaceColor]=deal('k');
errorbar([1,2,4,5],mm,[pci33b(1),pci36(1),pci63(1),pci66b(1)].*100-mm,'k.');
ylabel('FC rate (%)')
set(gca,'XTick',[1,2,4,5],'XTickLabel',{'3s duration & 3s sample',...
    '3s duration & 6s sample',...
    '6s duration & 3s sample',...
    '6s duration & 6s sample'},'XTickLabelRotation',45,'FontSize',10)

[~,~,p]=crosstab([ones(nnz_pair_33b,1); ...
    2*ones(nnz_pair_36,1)],...
    [(1:nnz_pair_33b)>nnz_sig_33b,...
    (1:nnz_pair_36)>nnz_sig_36]);
disp(p);

[~,~,p]=crosstab([ones(nnz_pair_63,1); ...
    2*ones(nnz_pair_66b,1)],...
    [(1:nnz_pair_63)>nnz_sig_63,...
    (1:nnz_pair_66b)>nnz_sig_66b]);
disp(p);
exportgraphics(fh,'duration_waveid_fc.pdf','ContentType','vector')


