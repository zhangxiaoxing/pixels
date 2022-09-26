function com_str_=get_pct_com_map(pct_meta,opt)
arguments
    pct_meta  (1,1) struct
    opt.onepath (1,:) char = '' % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
    opt.one_SU_showcase (1,1) logical = false % for the TCOM-FC joint showcase
    opt.append_late_delay (1,1) logical = false % Uses stats from early delay but include illustration for late delay
    opt.band_width (1,1) double {mustBeMember(opt.band_width,1:2)} = 1
end

%TODO proper declaration
persistent com_str opt_ pct_meta_

if isempty(com_str) || ~isequaln(opt,opt_) || ~isequaln(pct_meta_,pct_meta) 
    meta=ephys.util.load_meta('skip_stats',true);
    if strlength(opt.onepath)==0
        usess=unique(meta.sess);
    else
        dpath=regexp(opt.onepath,'(?<=SPKINFO[\\/]).*$','match','once');
        if isempty(dpath)
            dpath=opt.onepath;
        end
        usess=ephys.path2sessid(dpath);
    end

    homedir=ephys.util.getHomedir('type','raw');
    com_str=struct();

    pct_sel=struct();

    %%  mixed
    pct_sel.s1d3=pct_meta.wave_id==1;
    pct_sel.s1d6=pct_meta.wave_id==2;
    pct_sel.s2d3=pct_meta.wave_id==3;
    pct_sel.s2d6=pct_meta.wave_id==4;

    pct_sel.olf_s1=pct_meta.wave_id==5;
    pct_sel.olf_s2=pct_meta.wave_id==6;

    pct_sel.dur_d3=pct_meta.wave_id==7;
    pct_sel.dur_d6=pct_meta.wave_id==8;

    for sessid=reshape(usess,1,[])
        fpath=fullfile(homedir,ephys.sessid2path(sessid),'FR_All_ 250.hdf5');
        sesssel=meta.sess==sessid;
        if ~any(sesssel), continue;end
        fpath=replace(fpath,'\',filesep());
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        trl.cs1d3=find(trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
        trl.es1d3=find(trials(:,5)==4 & trials(:,8)==3 & trials(:,10)==0);
        trl.cs2d3=find(trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
        trl.es2d3=find(trials(:,5)==8 & trials(:,8)==3 & trials(:,10)==0);
        trl.cs1d6=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);
        trl.es1d6=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,10)==0);
        trl.cs2d6=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);
        trl.es2d6=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,10)==0);

        sess=['s',num2str(sessid)];


        for ff=["s1d3","s2d3","s1d6","s2d6","olf_s1","olf_s2","dur_d3","dur_d6"]
            mcid1=meta.allcid(sesssel & pct_sel.(ff));
            [~,msel1]=ismember(mcid1,suid);

            %         if opt.rnd_half
            %             for ff=["c1a","c1b","c1e"]
            %                 com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            %             end
            %             if opt.curve
            %                 for ff=["c1aheat","c1acurve","c1aanticurve",...
            %                         "c1bheat","c1bcurve","c1banticurve",...
            %                         "c1eheat","c1ecurve","c1eanticurve"]
            %                     com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            %                 end
            %             end
            %         else

            com_str.(['s',num2str(sessid)]).(ff).com=containers.Map('KeyType','int32','ValueType','any');
            com_str.(['s',num2str(sessid)]).(ff).com4plot=containers.Map('KeyType','int32','ValueType','any');

            if opt.curve
%                 if startsWith(ff,'s')
                for cc=["s1d3","s2d3","s1d6","s2d6"]
                    com_str.(['s',num2str(sessid)]).(ff).(cc)=containers.Map('KeyType','int32','ValueType','any');
                end
%                 elseif startsWith(ff,'o')
%                     for cc=["olf_s1","olf_s2"]
%                         com_str.(['s',num2str(sessid)]).(ff).(cc)=containers.Map('KeyType','int32','ValueType','any');
%                     end
%                 else
%                     for cc=["dur_d3","dur_d6"]
%                         com_str.(['s',num2str(sessid)]).(ff).(cc)=containers.Map('KeyType','int32','ValueType','any');
%                     end
%                 end
            end
            if isempty(msel1)
                continue
            end
            assert(~opt.rnd_half,'Random half-half not implemented yet')
            %         if opt.rnd_half
            %             c1a=randsample(c1sel,floor(numel(c1sel)./2));
            %             c1b=c1sel(~ismember(c1sel,c1a));
            %             if nnz(c1a)>2 && nnz(c1b)>2 && nnz(e1sel)>2
            %                 com_str=per_su_process(sess,suid,msel1,fr,trl,com_str,'c1a',opt);
            %                 com_str=per_su_process(sess,suid,msel1,fr,trl,com_str,'c1b',opt);
            %                 com_str=per_su_process(sess,suid,msel1,fr,trl,com_str,'c1e',opt);
            %             end
            %         else
            %         end
            if ismember(ff,["s1d3","s2d3","s1d6","s2d6"])
                com_str=per_su_process(sess,suid,msel1,fr,trl,com_str,ff,opt);
            elseif ismember(ff,["olf_s1","olf_s2"])
                com_str=per_su_process_olf(sess,suid,msel1,fr,trl,com_str,ff,opt);
            else 
                com_str=per_su_process_dur(sess,suid,msel1,fr,trl,com_str,ff,opt);
            end
        end
    end
end
com_str_=com_str;
opt_=opt;
pct_meta_=pct_meta;
end

function com_str=per_su_process(sess,suid,msel,fr,trls,com_str,type,opt)
    for su=reshape(msel,1,[])
        classmm=[];
        for ff=["cs1d3","cs2d3","cs1d6","cs2d6"]
            ffmat=squeeze(fr(trls.(ff),su,:));
            classmm=cat(1,classmm,mean(ffmat(:,17:40),1));
        end
        basemm=mean([classmm(:,1:12)],'all');
        S=max(([classmm(:,1:12)]-basemm),[],'all');% removed abs

        classnn=(classmm-basemm)./S;
        [~,typeidx]=ismember("c"+type,["cs1d3","cs2d3","cs1d6","cs2d6"]);
        mm_pref=classnn(typeidx,1:12);
        mm_pref(mm_pref<0)=0;
        if ~any(mm_pref>0)
            disp(strjoin({sess,num2str(suid(su)),char(type),'PEAK mismatch, TCOM set to -1'},','))
            keyboard();
            continue
        else
            com=sum((1:12).*mm_pref)./sum(mm_pref);
            com_str.(sess).(type).com(suid(su))=com;
        end

        if opt.curve
            basemm=mean([classmm(:,1:12);classmm(3:4,13:24)],'all');
            S=max(([classmm(:,1:12);classmm(3:4,13:24)]-basemm),[],'all');% removed abs
            classnn=(classmm-basemm)./S;
            for typeidx=1:4
                cc=subsref(["s1d3","s2d3","s1d6","s2d6"],struct(type='()',subs={{typeidx}}));
                if typeidx<3
                    com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,1:12);
                else
                    com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,:);
                end
            end
            if contains(type,'d3')
                com_str.(sess).(type).com4plot(suid(su))=com;
            else
                if contains(type,'s1')
                    ppref=classnn(3,:);
                else
                    ppref=classnn(4,:);
                end
                ppref(ppref<0)=0;
                if any(ppref>0)
                    com=sum((1:24).*ppref)./sum(ppref);
                    com_str.(sess).(type).com4plot(suid(su))=com;
                end
            end
        end
    end
end


function com_str=per_su_process_olf(sess,suid,msel,fr,trls,com_str,type,opt)

    for su=reshape(msel,1,[])
        classmm=[];
        for ff=["cs1d3","cs2d3","cs1d6","cs2d6"]
            ffmat=squeeze(fr(trls.(ff),su,:));
            classmm=cat(1,classmm,mean(ffmat(:,17:40),1));
        end
        basemm=mean([classmm(:,1:12)],'all');
        S=max(([classmm(:,1:12)]-basemm),[],'all');% removed abs
        classnn=(classmm-basemm)./S;
        
        if contains(type,'s1')
            mm_pref=mean(classnn([1,3],1:12));
        else
            mm_pref=mean(classnn([2 4],1:12));
        end
        mm_pref(mm_pref<0)=0;
        if ~any(mm_pref>0)
            keyboard();
            continue
        else
            com=sum((1:12).*mm_pref)./sum(mm_pref);
            com_str.(sess).(type).com(suid(su))=com;
        end
        
        if opt.curve

            basemm=mean([classmm(:,1:12);classmm(3:4,13:24)],'all');
            S=max(([classmm(:,1:12);classmm(3:4,13:24)]-basemm),[],'all');% removed abs
            classnn=(classmm-basemm)./S;

            for typeidx=1:4
                cc=subsref(["s1d3","s2d3","s1d6","s2d6"],struct(type='()',subs={{typeidx}}));
                if typeidx<3
                    com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,1:12);
                else
                    com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,:);
                end
            end

            if contains(type,'s1')
                ppref=classnn(3,:);
            else
                ppref=classnn(4,:);
            end
            ppref(ppref<0)=0;
            if any(ppref>0)
                com=sum((1:24).*ppref)./sum(ppref);
                com_str.(sess).(type).com4plot(suid(su))=com;
            end
        end
    end
end

function com_str=per_su_process_dur(sess,suid,msel,fr,trls,com_str,type,opt)
    for su=reshape(msel,1,[])
        classmm=[];
        for ff=["cs1d3","cs2d3","cs1d6","cs2d6"]
            ffmat=squeeze(fr(trls.(ff),su,:));
            classmm=cat(1,classmm,mean(ffmat(:,17:40),1));
        end
        basemm=mean([classmm(:,1:12)],'all');
        S=max(([classmm(:,1:12)]-basemm),[],'all');% removed abs
        classnn=(classmm-basemm)./S;

        if contains(type,'d3')
            mm_pref=mean(classnn(1:2,1:12));
        else
            mm_pref=mean(classnn(3:4,1:12));
        end
        mm_pref(mm_pref<0)=0;
        if ~any(mm_pref>0)
            disp(strjoin({sess,num2str(suid(su)),char(type),'PEAK mismatch, TCOM set to -1'},','))
            keyboard();
            continue
        else
            com=sum((1:12).*mm_pref)./sum(mm_pref);
            com_str.(sess).(type).com(suid(su))=com;
        end
        
        if opt.curve
            basemm=mean([classmm(:,1:12);classmm(3:4,13:24)],'all');
            S=max(([classmm(:,1:12);classmm(3:4,13:24)]-basemm),[],'all');% removed abs
            classnn=(classmm-basemm)./S;
            for typeidx=1:4
                cc=subsref(["s1d3","s2d3","s1d6","s2d6"],struct(type='()',subs={{typeidx}}));
                if typeidx<3
                    com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,1:12);
                else
                    com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,:);
                end
            end

            if contains(type,'d3')
                ppref=mean(classnn(1:2,1:12));
                ppref(ppref<0)=0;
                if any(ppref>0)
                    com=sum((1:12).*ppref)./sum(ppref);
                    com_str.(sess).(type).com4plot(suid(su))=com;
                end
            else
                ppref=mean(classnn(3:4,:));
                ppref(ppref<0)=0;
                if any(ppref>0)
                    com=sum((1:24).*ppref)./sum(ppref);
                    com_str.(sess).(type).com4plot(suid(su))=com;
                end
            end

        end
    end
end

