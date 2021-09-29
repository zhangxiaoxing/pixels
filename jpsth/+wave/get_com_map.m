%TODO Error trial
function com_str_=get_com_map(opt)
arguments
    opt.onepath (1,:) char = '' % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.per_sec_stats (1,1) logical = false % calculate COM using per-second mean as basis for normalized firing rate, default is coss-delay mean
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
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
        fpath=replace(fpath,'\',filesep());
        pc_stem=replace(dpath,'/','\');
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
        s1sel=find(trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0);
        s2sel=find(trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0);
        
        e1sel=find(trial(:,5)==4 & trial(:,8)==6 & trial(:,10)==0);
        e2sel=find(trial(:,5)==8 & trial(:,8)==6 & trial(:,10)==0);
        
        sess=['s',num2str(sessid)];
        %     if sum(trial(:,9))<40,continue;end %meta data obtained from processed
        %     welltrained dataset
        if opt.rnd_half
            for ff=["s1a","s2a","s1b","s2b","s1e","s2e"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["s1aheat","s2aheat","s1acurve","s2acurve",...
                        "s1bheat","s2bheat","s1bcurve","s2bcurve",...
                        "e1heat","e2heat","e1curve","e2curve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
            s1a=randsample(s1sel,floor(numel(s1sel)./2));
            s1b=s1sel(~ismember(s1sel,s1a));
            s2a=randsample(s2sel,floor(numel(s2sel)./2));
            s2b=s2sel(~ismember(s2sel,s2a));
            if nnz(s1a)>2 && nnz(s1b)>2 && nnz(s2a)>2 && nnz(s2b)>2 && nnz(e1sel)>2 && nnz(e2sel)>2
                com_str=per_su_process(sess,suid,msel1,fr,s1a,s2a,com_str,'s1a',opt);
                com_str=per_su_process(sess,suid,msel1,fr,s1b,s2b,com_str,'s1b',opt);
                com_str=per_su_process(sess,suid,msel2,fr,s2a,s1a,com_str,'s2a',opt);
                com_str=per_su_process(sess,suid,msel2,fr,s2b,s1b,com_str,'s2b',opt);
                com_str=per_su_process(sess,suid,msel1,fr,e1sel,e2sel,com_str,'s1e',opt);
                com_str=per_su_process(sess,suid,msel2,fr,e2sel,e1sel,com_str,'s2e',opt);
            end
        else
            for ff=["s1","s2"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["s1heat","s2heat","s1curve","s2curve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
            
            com_str=per_su_process(sess,suid,msel1,fr,s1sel,s2sel,com_str,'s1',opt);
            com_str=per_su_process(sess,suid,msel2,fr,s2sel,s1sel,com_str,'s2',opt);
        end
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
    perfmat=squeeze(fr(pref_sel,su,stats_window));
    npmat=squeeze(fr(nonpref_sel,su,stats_window));
    basemm=mean([mean(perfmat,1);mean(npmat,1)]);
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
    if opt.curve
        curve=mm_pref;
        com_str.(sess).([samp,'curve'])(suid(su))=curve;
        heatcent=squeeze(fr(pref_sel,su,stats_window))-basemm; %centralized norm. firing rate for heatmap plot
        heatnorm=heatcent./max(heatcent);
        if opt.rnd_half
            heatnorm=mean(heatnorm);
        else
            heatnorm(heatnorm<0)=0;
            if size(heatnorm,1)>10
                cc=arrayfun(@(x) min(corrcoef(heatnorm(x,:),curve),[],'all'),1:size(heatnorm,1));
                [~,idx]=sort(cc,'descend');
                heatnorm=heatnorm(idx(1:10),:);
            end
        end
        com_str.(sess).([samp,'heat'])(suid(su))=heatnorm;
    end
end
end
