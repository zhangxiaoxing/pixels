function com_str_=get_com_map(opt)
arguments
    opt.onepath (1,:) char = '' % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.per_sec_stats (1,1) logical = false % calculate COM using per-second mean as basis for normalized firing rate, default is coss-delay mean
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
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
%TODO if decision
if opt.decision
    stats_window=41:44;
else
    stats_window=17:40;
end
for su=reshape(msel,1,[])
%     if strcmp(sess,'s67') && suid(su)==10993,keyboard();end
    basemm=mean([mean(squeeze(fr(pref_sel,su,stats_window)));mean(squeeze(fr(nonpref_sel,su,stats_window)))]);
    if ~opt.per_sec_stats
        basemm=mean(basemm);
    end
    %TODO check the effect of smooth
    mm=smooth(squeeze(mean(fr(pref_sel,su,:))),5).';
    mm_pref=mm(stats_window)-basemm;
    if max(mm_pref)<=0,continue;end
    mm_pref=mm_pref./max(mm_pref);
    mm_pref(mm_pref<0)=0;
    com=sum((1:numel(stats_window)).*mm_pref)./sum(mm_pref);
    com_str.(sess).(samp)(suid(su))=com;
    curve=mm_pref;
    com_str.(sess).([samp,'curve'])(suid(su))=curve;
        
    heatcent=squeeze(fr(pref_sel,su,stats_window))-basemm; %centralized norm. firing rate for heatmap plot
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