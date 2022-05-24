%TODO to reuse wave.get_comm_map and purge duplicated codes

function [fh,data]=plot_wave(opt)
arguments
    opt.plot_error (1,1) logical = false
    opt.per_sec_stats (1,1) logical = true
    opt.mtype (1,:) char {mustBeMember(opt.mtype,{'memory','nonmem'})}='memory'
end
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
    
    %TODO nonmemory, incongruent?
    mcid1=meta_str.allcid(meta_str.mem_type==2 & sesssel.');
    msel1=find(ismember(suid,mcid1));
    mcid2=meta_str.allcid(meta_str.mem_type==4 & sesssel.');
    msel2=find(ismember(suid,mcid2));

    s1sel=trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
    s2sel=trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;

    if opt.plot_error
        es1sel=trial(:,5)==4 & trial(:,8)==6 & trial(:,10)==0;
        es2sel=trial(:,5)==8 & trial(:,8)==6 & trial(:,10)==0;
        if nnz(es1sel)<2 || nnz(es2sel)<2
            continue
        end
    end
    
    for su=reshape(msel1,1,[])
        basemm=mean([mean(squeeze(fr(s1sel,su,:)));mean(squeeze(fr(s2sel,su,:)))]);
        basestd=mean([std(squeeze(fr(s1sel,su,:)));std(squeeze(fr(s2sel,su,:)))]);
        basestd(basestd==0 & basemm==0)=1;
        if ~opt.per_sec_stats
            basemm=mean(basemm(17:40));
            basestd=mean(basestd(17:40));
        end
        if any(basestd==0),keyboard();end
        if opt.plot_error
            h11=[h11;((squeeze(mean(fr(es1sel,su,:))).'-basemm)./basestd)];
            h21=[h21;((squeeze(mean(fr(es2sel,su,:))).'-basemm)./basestd)];
        else
            h11=[h11;((squeeze(mean(fr(s1sel,su,:))).'-basemm)./basestd)];   
            h21=[h21;((squeeze(mean(fr(s2sel,su,:))).'-basemm)./basestd)];  
        end
    end
    
    for su=reshape(msel2,1,[])
        basemm=mean([mean(squeeze(fr(s1sel,su,:)));mean(squeeze(fr(s2sel,su,:)))]);
        basestd=mean([std(squeeze(fr(s1sel,su,:)));std(squeeze(fr(s2sel,su,:)))]);
        basestd(basestd==0 & basemm==0)=1;
        if any(basestd==0),keyboard();end
        if ~opt.per_sec_stats
            basemm=mean(basemm(17:40));
            basestd=mean(basestd(17:40));
        end
        if opt.plot_error
            h12=[h12;((squeeze(mean(fr(es1sel,su,:))).'-basemm)./basestd)];
            h22=[h22;((squeeze(mean(fr(es2sel,su,:))).'-basemm)./basestd)];
        else
            h12=[h12;((squeeze(mean(fr(s1sel,su,:))).'-basemm)./basestd)];
            h22=[h22;((squeeze(mean(fr(s2sel,su,:))).'-basemm)./basestd)];
        end
    end
end

h11pos=h11;h11pos(h11pos<0)=0;
h22pos=h22;h22pos(h22pos<0)=0;
com1=sum(((1:24).*h11pos(:,17:40)),2)./sum(h11pos(:,17:40),2);
com2=sum(((1:24).*h22pos(:,17:40)),2)./sum(h22pos(:,17:40),2);

[scom1,iidx1]=sort(com1);
[scom2,iidx2]=sort(com2);

fh=figure('Color','w','Position',[100,100,500,600]);
wave.plotOne(h11(iidx1,:),1,'S1 trials','com',scom1+16);
wave.plotOne(h21(iidx1,:),2,'S2 trials');
wave.plotOne(h12(iidx2,:),3,'S1 trials');
wave.plotOne(h22(iidx2,:),4,'S2 trials','com',scom2+16);
[data.h11,data.h21,data.h12,data.h22]=deal(h11,h21,h12,h22);
end