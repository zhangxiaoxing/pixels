function com_str_=get_com_map(opt)
arguments
    opt.onepath (1,:) char = ''
    opt.peak (1,1) logical = false
    opt.curve (1,1) logical = false
    opt.per_sec_stats (1,1) logical = false
    opt.selidx_curve (1,1) logical = false
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
        %TODO nonmem,incongruent
        mcid1=meta_str.allcid(meta_str.mem_type==2 & sesssel.');
        msel1=find(ismember(suid,mcid1));
        mcid2=meta_str.allcid(meta_str.mem_type==4 & sesssel.');
        msel2=find(ismember(suid,mcid2));
        if isempty(msel1) && isempty(msel2), continue; end
        sessid=ephys.path2sessid(pc_stem);
        com_str.(['s',num2str(sessid)]).s1=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(sessid)]).s2=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(sessid)]).s1heat=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(sessid)]).s2heat=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(sessid)]).s1curve=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(sessid)]).s2curve=containers.Map('KeyType','int32','ValueType','any');
        %     if sum(trial(:,9))<40,continue;end %meta data obtained from processed
        %     welltrained dataset
        s1sel=trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
        s2sel=trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
        sess=['s',num2str(sessid)];
        com_str=per_su_process(sess,suid,msel1,fr,s1sel,s2sel,com_str,'s1',opt);
        com_str=per_su_process(sess,suid,msel2,fr,s2sel,s1sel,com_str,'s2',opt);

        if ~strlength(opt.onepath)==0
            break;
        end
    end
end
com_str_=com_str;
end

function com_str=per_su_process(sess,suid,msel,fr,pref_sel,nonpref_sel,com_str,samp,opt)
for su=reshape(msel,1,[])
    basemm=mean([mean(squeeze(fr(pref_sel,su,17:40)));mean(squeeze(fr(nonpref_sel,su,17:40)))]);
    if ~opt.per_sec_stats
        basemm=mean(basemm);
    end
    mm=smooth(squeeze(mean(fr(pref_sel,su,:))),5).';
    mm_pref=mm(17:40)-basemm;
    mm_pref=mm_pref./max(mm_pref);
    mm_pref(mm_pref<0)=0;
    if opt.peak
        [~,pidx]=max(mm_pref);
        com_str.(sess).(samp)(suid(su))=pidx;
    else
        com=sum((1:24).*mm_pref)./sum(mm_pref);
        com_str.(sess).(samp)(suid(su))=com;
    end
    
    if opt.selidx_curve
        fr_pref=squeeze(mean(fr(pref_sel,su,17:40)));
        fr_nonp=squeeze(mean(fr(nonpref_sel,su,17:40)));
        curve=(fr_pref-fr_nonp)./(fr_pref+fr_nonp);
        com_str.(sess).([samp,'curve'])(suid(su))=curve;
            
    else
        curve=mm_pref;
        com_str.(sess).([samp,'curve'])(suid(su))=curve;
            
    end
    
    heatcent=squeeze(fr(pref_sel,su,17:40))-basemm;
    heatcent(heatcent<0)=0;
    heatnorm=heatcent./max(heatcent);
    if size(heatnorm,1)>10
        cc=arrayfun(@(x) min(corrcoef(heatnorm(x,:),curve),[],'all'),1:size(heatnorm,1));
        [~,idx]=sort(cc,'descend');
        heatnorm=heatnorm(idx(1:10),:);
    end
    com_str.(sess).([samp,'heat'])(suid(su))=heatnorm;    
end
end