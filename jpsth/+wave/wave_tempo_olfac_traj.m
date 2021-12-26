% FR->D3S1 D6S1 D3S2 D6S2
[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir('type','raw');
pca_str=struct();
pcamat=[];
for ii=reshape(cell2mat(sessmap.keys()),1,[])
    disp(ii)
    fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
    fr=h5read(fpath,'/FR_All');
    trial=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');
    d3s1sel=find(trial(:,5)==4 & trial(:,8)==3 & trial(:,9)>0 & trial(:,10)>0);
    d3s2sel=find(trial(:,5)==8 & trial(:,8)==3 & trial(:,9)>0 & trial(:,10)>0);
    d6s1sel=find(trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0);
    d6s2sel=find(trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0);

%     e1sel=find(trial(:,5)==4 & trial(:,8)==opt.delay & trial(:,10)==0);
%     e2sel=find(trial(:,5)==8 & trial(:,8)==opt.delay & trial(:,10)==0);

    sess=['s',num2str(ii)];
    pca_str.(['s',num2str(ii)])=struct();
    for ff=["d3s1","d3s2","d6s1","d6s2"]
        pca_str.(['s',num2str(ii)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
    end
    pca_str=per_su_process(sess,suid,fr,d3s1sel,d3s2sel,d6s1sel,d6s2sel,pca_str);
    % totally unnecessary detour
    pcamat=[pcamat;cell2mat(pca_str.(['s',num2str(ii)]).d3s1.values(num2cell(suid.'))).',...
        cell2mat(pca_str.(['s',num2str(ii)]).d3s2.values(num2cell(suid.'))).',...
        cell2mat(pca_str.(['s',num2str(ii)]).d6s1.values(num2cell(suid.'))).',...
        cell2mat(pca_str.(['s',num2str(ii)]).d6s2.values(num2cell(suid.'))).'];
end

% magicnum=56;
normmat=normalize(pcamat.','range');
%PCA
[coef,score,latent]=pca(normmat);
fh=figure('Color','w');
hold on
d3s1h=plot3(score(1:56,1),score(1:56,2),score(1:56,3),'--r');
d3s2h=plot3(score((1:56)+56,1),score((1:56)+56,2),score((1:56)+56,3),'--b');
d6s1h=plot3(score((1:56)+112,1),score((1:56)+112,2),score((1:56)+112,3),'-r');
d6s2h=plot3(score((1:56)+168,1),score((1:56)+168,2),score((1:56)+168,3),'-b');
for tt=2:2:14
    text(score(tt*4,1),score(tt*4,2),score(tt*4,3)-4,num2str(tt-4),'FontSize',12);
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
legend([d6s1h,d6s2h,d3s1h,d3s2h],{'S1 6s','S2 6s','S1 3s','S2 3s'},Location='northoutside',Orientation='horizontal')
v=VideoWriter('tempo_olfac_pca.mp4','MPEG-4');
open(v)
for az=20:0.5:380
    view(az,10)
    writeVideo(v,getframe(fh))
end
close(v)


% coding direction
% odor cd
ocd=mean(normmat([17:28,(17:40)+112],:))-mean(normmat([(17:28)+56,(17:40)+168],:));
ocd=ocd/norm(ocd);
% duration cd
dcd=mean(normmat([1:12,(1:12)+56],:))-mean(normmat([(1:12)+112,(1:12)+168],:));
dcd=dcd/norm(dcd);
fh=figure('Color','w');
hold on
grid on
d3s1h=plot3(normmat(1:56,:)*(ocd.'),normmat(1:56,:)*(dcd.'),1:56,'--r');
d3s2h=plot3(normmat((1:56)+56,:)*(ocd.'),normmat((1:56)+56,:)*(dcd.'),1:56,'--b');
d6s1h=plot3(normmat((1:56)+112,:)*(ocd.'),normmat((1:56)+112,:)*(dcd.'),1:56,'-r');
d6s2h=plot3(normmat((1:56)+168,:)*(ocd.'),normmat((1:56)+168,:)*(dcd.'),1:56,'-b');

text(normmat(20,:)*(ocd.'),normmat(20,:)*(dcd.'),20,'S1,3s','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle');
text(normmat(24+56,:)*(ocd.'),normmat(24+56,:)*(dcd.'),24,'S2,3s','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle');
text(normmat(28+112,:)*(ocd.'),normmat(28+112,:)*(dcd.'),28,'S1,6s','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle');
text(normmat(32+168,:)*(ocd.'),normmat(32+168,:)*(dcd.'),32,'S2,6s','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle');

set(gca,'ZTick',8:8:56,'ZTickLabel',-2:2:10)
xlabel('Olfactory CD')
ylabel('Duration prediction CD')
zlabel('Time (s)')
% legend([d6s1h,d6s2h,d3s1h,d3s2h],{'S1 6s','S2 6s','S1 3s','S2 3s'},Location='northoutside',Orientation='horizontal')
v=VideoWriter('tempo_olfac_cd.mp4','MPEG-4');
open(v)
for az=-30:0.5:330
    view(az,15)
    writeVideo(v,getframe(fh))
end
close(v)




function pca_str=per_su_process(sess,suid,fr,d3s1sel,d3s2sel,d6s1sel,d6s2sel,pca_str)
% arguments
%     opt.msel double =[]
% end

for su=1:size(fr,2)
%     d3s1fr=squeeze(fr(d3s1sel,su,:));
%     d3s2fr=squeeze(fr(d3s2sel,su,:));
%     d6s1fr=squeeze(fr(d6s1sel,su,:));
%     d6s2fr=squeeze(fr(d6s2sel,su,:));

    pca_str.(sess).d3s1(suid(su))=squeeze(mean((fr(d3s1sel,su,:))));
    pca_str.(sess).d3s2(suid(su))=squeeze(mean((fr(d3s2sel,su,:))));
    pca_str.(sess).d6s1(suid(su))=squeeze(mean((fr(d6s1sel,su,:))));
    pca_str.(sess).d6s2(suid(su))=squeeze(mean((fr(d6s2sel,su,:))));
end
end


