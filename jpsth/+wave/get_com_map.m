function com_str_=get_com_map(opt)
arguments
    opt.onepath (1,:) char = ''
end
persistent com_str onepath_
if isempty(onepath_), onepath_='';end
if isempty(com_str) || ~strcmp(opt.onepath, onepath_)
    meta_str=ephys.util.load_meta('type','neupix');
    homedir=ephys.util.getHomedir('type','raw');
    fl=dir(fullfile(homedir,'**','FR_All_ 250.hdf5'));
    com_str=struct();
    for ii=1:size(fl,1)
        if strcmp(opt.onepath,'')
            dpath=regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once');
        else
            dpath=opt.onepath;
        end
        pc_stem=replace(dpath,'/',filesep());
        sesssel=startsWith(meta_str.allpath,pc_stem);
        if ~any(sesssel), continue;end
        fr=h5read(fullfile(fl(ii).folder,fl(ii).name),'/FR_All');
        trial=h5read(fullfile(fl(ii).folder,fl(ii).name),'/Trials');
        suid=h5read(fullfile(fl(ii).folder,fl(ii).name),'/SU_id');
        mcid1=meta_str.allcid(meta_str.mem_type==2 & sesssel.');
        msel1=find(ismember(suid,mcid1));
        mcid2=meta_str.allcid(meta_str.mem_type==4 & sesssel.');
        msel2=find(ismember(suid,mcid2));
        if isempty(msel1) && isempty(msel2), continue; end
        sessid=ephys.path2sessid(pc_stem);
        com_str.(['s',num2str(sessid)])=containers.Map('KeyType','int32','ValueType','any');
        %     if sum(trial(:,9))<40,continue;end %meta data obtained from processed
        %     welltrained dataset
        s1sel=trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
        s2sel=trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
        
        for su=reshape(msel1,1,[])
            mm1=squeeze(fr(s1sel,su,17:40));
%             mm2=squeeze(fr(s2sel,su,17:40));
%             mm=mean([mean(mm1);mean(mm2)]);
%             stdd=mean([std(mm1);std(mm2)]);
%             norm=(mean(mm1)-mm)./stdd;
%             norm=mean(mm1)-mean(mm2);
%             norm=norm-min(norm);
            norm=mean(mm1);
            com=sum((1:24).*norm)./sum(norm);
            com_str.(['s',num2str(sessid)])(suid(su))=com;
        end
        
        for su=reshape(msel2,1,[])
%             mm1=squeeze(fr(s1sel,su,17:40));
            mm2=squeeze(fr(s2sel,su,17:40));
%             mm=mean([mean(mm1);mean(mm2)]);
%             stdd=mean([std(mm1);std(mm2)]);
%             norm=(mean(mm2)-mm)./stdd;
%             norm=mean(mm2)-mean(mm1);
%             norm=norm-min(norm);
            norm=mean(mm2);
            com=sum((1:24).*norm)./sum(norm);
            com_str.(['s',num2str(sessid)])(suid(su))=com;
        end
        if ~strcmp(opt.onepath,'')
            break;
        end
    end
end
com_str_=com_str;
end