meta_str=ephys.util.load_meta('type','neupix');
homedir=ephys.util.getHomedir('type','raw');
fl=dir(fullfile(homedir,'**','FR_All_ 250.hdf5'));
h1=[];
h2=[];
for ii=1:size(fl,1)
    pc_stem=replace(regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once'),'/','\');
    sesssel=startsWith(meta_str.allpath,pc_stem);
    if ~any(sesssel), continue;end
    fr=h5read(fullfile(fl(ii).folder,fl(ii).name),'/FR_All');
    trial=h5read(fullfile(fl(ii).folder,fl(ii).name),'/Trials');
    suid=h5read(fullfile(fl(ii).folder,fl(ii).name),'/SU_id');
    mcid=meta_str.allcid(meta_str.mem_type>2 & sesssel.');
    msel=ismember(suid,mcid);
%     sel_id=
%     if sum(trial(:,9))<40,continue;end
    s1sel=trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
    s2sel=trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
    h1=[h1;reshape(squeeze(mean(fr(s1sel,msel,:))),[],size(fr,3))];
    h2=[h2;reshape(squeeze(mean(fr(s2sel,msel,:))),[],size(fr,3))];
end

normh1=h1./max(h1,[],2);
normh2=h2./max(h2,[],2);
std1=std(normh1,0,2);
std2=std(normh2,0,2);

normh1=normh1(std1>0 & std2>0,:);
normh2=normh2(std1>0 & std2>0,:);

com1=sum(((1:24).*normh1(:,17:40)),2)./sum(normh1(:,17:40),2);
com2=sum(((1:24).*normh2(:,17:40)),2)./sum(normh2(:,17:40),2);

[~,iidx1]=sort(com1);
[~,iidx2]=sort(com2);

figure('Color','w');
subplot(1,2,1)
imagesc(normh1(iidx1,:))
colormap('jet')
ax=gca();
ax.YDir='normal';

subplot(1,2,2)
imagesc(normh2(iidx1,:))
colormap('jet')
ax=gca();
ax.YDir='normal';

figure('Color','w');
subplot(1,2,1)
imagesc(normh1(iidx2,:))
colormap('jet')
ax=gca();
ax.YDir='normal';

subplot(1,2,2)
imagesc(normh2(iidx2,:))
colormap('jet')
ax=gca();
ax.YDir='normal';

