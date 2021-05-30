meta_str=ephys.util.load_meta('type','neupix');
homedir=ephys.util.getHomedir('type','raw');
fl=dir(fullfile(homedir,'**','FR_All_ 250.hdf5'));
h11=[];
h21=[];
h12=[];
h22=[];

for ii=1:size(fl,1)
    pc_stem=replace(regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once'),'/','\');
    sesssel=startsWith(meta_str.allpath,pc_stem);
    if ~any(sesssel), continue;end
    fr=h5read(fullfile(fl(ii).folder,fl(ii).name),'/FR_All');
    trial=h5read(fullfile(fl(ii).folder,fl(ii).name),'/Trials');
    suid=h5read(fullfile(fl(ii).folder,fl(ii).name),'/SU_id');
    mcid1=meta_str.allcid(meta_str.mem_type==2 & sesssel.');
    msel1=find(ismember(suid,mcid1));
    mcid2=meta_str.allcid(meta_str.mem_type==4 & sesssel.');
    msel2=find(ismember(suid,mcid2));
%     sel_id=
%     if sum(trial(:,9))<40,continue;end
    s1sel=trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
    s2sel=trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;

    for su=reshape(msel1,1,[])
        basemm=mean(fr(s1sel | s2sel,su,17:40),'all');
        basestd=std(fr(s1sel | s2sel,su,17:40),0,'all');
        h11=[h11;((squeeze(mean(fr(s1sel,su,:)))-basemm)./basestd).'];    
        h21=[h21;((squeeze(mean(fr(s2sel,su,:)))-basemm)./basestd).'];    
    end
    
    for su=reshape(msel2,1,[])
        basemm=mean(fr(s1sel | s2sel,su,17:40),'all');
        basestd=std(fr(s1sel | s2sel,su,17:40),0,'all');
        h12=[h12;((squeeze(mean(fr(s1sel,su,:)))-basemm)./basestd).'];    
        h22=[h22;((squeeze(mean(fr(s2sel,su,:)))-basemm)./basestd).'];    
    end
end

h11pos=h11-min(h11,[],'all');
h22pos=h22-min(h22,[],'all');
com1=sum(((1:24).*h11pos(:,17:40)),2)./sum(h11pos(:,17:40),2);
com2=sum(((1:24).*h22pos(:,17:40)),2)./sum(h22pos(:,17:40),2);

[~,iidx1]=sort(com1);
[~,iidx2]=sort(com2);

fh=figure('Color','w');
wave.plotOne(h11(iidx1,:),1,'S1 trials');
wave.plotOne(h21(iidx1,:),2,'S2 trials');
wave.plotOne(h12(iidx2,:),3,'S1 trials');
wave.plotOne(h22(iidx2,:),4,'S2 trials');
