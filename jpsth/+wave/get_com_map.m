function com_str_=get_com_map(opt)
arguments
    opt.onepath (1,:) char = ''
    opt.peak (1,1) logical = false
    opt.curve (1,1) logical = false
    opt.per_sec_stats (1,1) logical = false
end
persistent com_str onepath_
if isempty(onepath_), onepath_='';end
if isempty(com_str) || ~strcmp(opt.onepath, onepath_)
    meta_str=ephys.util.load_meta('type','neupix');
    homedir=ephys.util.getHomedir('type','raw');
    fl=dir(fullfile(homedir,'**','FR_All_ 250.hdf5'));
    com_str=struct();
    for ii=1:size(fl,1)
        if strlength(opt.onepath)==0
            dpath=regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once');
            fpath=fullfile(fl(ii).folder,fl(ii).name);
        else
            dpath=regexp(opt.onepath,'(?<=SPKINFO[\\/]).*$','match','once');
            fpath=fullfile(homedir,dpath,'FR_All_ 250.hdf5');
        end
        pc_stem=replace(dpath,'/',filesep());
        sesssel=startsWith(meta_str.allpath,pc_stem);
        if ~any(sesssel), continue;end
        fr=h5read(fpath,'/FR_All');
        trial=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        mcid1=meta_str.allcid(meta_str.mem_type==2 & sesssel.');
        msel1=find(ismember(suid,mcid1));
        mcid2=meta_str.allcid(meta_str.mem_type==4 & sesssel.');
        msel2=find(ismember(suid,mcid2));
        if isempty(msel1) && isempty(msel2), continue; end
        sessid=ephys.path2sessid(pc_stem);
        com_str.(['s',num2str(sessid)]).s1=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(sessid)]).s2=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(sessid)]).s1curve=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(sessid)]).s2curve=containers.Map('KeyType','int32','ValueType','any');
        
        
        %     if sum(trial(:,9))<40,continue;end %meta data obtained from processed
        %     welltrained dataset
        s1sel=trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
        s2sel=trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
        

        for su=reshape(msel1,1,[])
            basemm=mean([mean(squeeze(fr(s1sel,su,:)));mean(squeeze(fr(s2sel,su,:)))]);
            basestd=mean([std(squeeze(fr(s1sel,su,:)));std(squeeze(fr(s2sel,su,:)))]);
            basestd(basestd==0 & basemm==0)=1;
            if ~opt.per_sec_stats
                basemm=mean(basemm(17:40));
                basestd=mean(basestd(17:40));
            else
                basemm=basemm(17:40);
                basestd=basestd(17:40);
            end
            mm1=((squeeze(mean(fr(s1sel,su,17:40)))-basemm)./basestd).';

            mm1(mm1<0)=0;
            if opt.peak
                [~,pidx]=max(mm1);
                com_str.(['s',num2str(sessid)]).s1(suid(su))=pidx;
            else
                com=sum((1:24).*mm1)./sum(mm1);
                com_str.(['s',num2str(sessid)]).s1(suid(su))=com;
            end
%             fr1=squeeze(mean(fr(s1sel,su,17:40)));
%             fr2=squeeze(mean(fr(s2sel,su,17:40)));
            com_str.(['s',num2str(sessid)]).s1curve(suid(su))=mm1;%(fr1-fr2)./(fr1+fr2);
        end
        for su=reshape(msel2,1,[])
            basemm=mean([mean(squeeze(fr(s1sel,su,:)));mean(squeeze(fr(s2sel,su,:)))]);
            basestd=mean([std(squeeze(fr(s1sel,su,:)));std(squeeze(fr(s2sel,su,:)))]);
            basestd(basestd==0 & basemm==0)=1;
            if ~opt.per_sec_stats
                basemm=mean(basemm(17:40));
                basestd=mean(basestd(17:40));
            else
                basemm=basemm(17:40);
                basestd=basestd(17:40);
            end
            mm2=((squeeze(mean(fr(s2sel,su,17:40)))-basemm)./basestd).';
            mm2(mm2<0)=0;
            if opt.peak
                [~,pidx]=max(mm2);
                com_str.(['s',num2str(sessid)]).s2(suid(su))=pidx;
            else
                com=sum((1:24).*mm2)./sum(mm2);
                com_str.(['s',num2str(sessid)]).s2(suid(su))=com;
            end
%             fr1=squeeze(mean(fr(s1sel,su,17:40)));
%             fr2=squeeze(mean(fr(s2sel,su,17:40)));
            com_str.(['s',num2str(sessid)]).s2curve(suid(su))=mm2;%(fr2-fr1)./(fr1+fr2);
            
        end
        if ~strlength(opt.onepath)==0
            break;
        end
    end
end
com_str_=com_str;
end