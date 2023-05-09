function [com_str_,com_str_h2]=get_pct_com_map(pct_meta,opt)
arguments
    pct_meta  (1,1) struct
    opt.one_sess (1,:) double = [] % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
    opt.err (1,1) logical = false %u stats in error trials
    opt.one_SU_showcase (1,1) logical = false % for the TCOM-FC joint showcase
    opt.append_late_delay (1,1) logical = false % Uses stats from early delay but include illustration for late delay
    opt.band_width (1,1) double {mustBeMember(opt.band_width,1:2)} = 1
    opt.early_smooth (1,1) logical = false % default value has changed as of 2023.03.01
end

%TODO proper declaration
persistent com_str opt_ pct_meta_

assert(~(opt.rnd_half && opt.err), 'Unable to handle 2Fold CV and error trial in same run')
if opt.rnd_half || isempty(com_str) || ~isequaln(opt,opt_) || ~isequaln(pct_meta_,pct_meta) 
    meta=ephys.util.load_meta('skip_stats',true);
    if isempty(opt.one_sess)
        usess=unique(meta.sess);
    else
        usess=opt.one_sess;
    end

    homedir=ephys.util.getHomedir('type','raw');
    com_str=struct();
    if opt.rnd_half
        com_str_h2=struct();
    end

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
%         if sessid==11
%             keyboard()
%         end
        fpath=fullfile(homedir,ephys.sessid2path(sessid),'FR_All_ 250.hdf5');
        sesssel=meta.sess==sessid;
        if ~any(sesssel), continue;end
        fpath=replace(fpath,'\',filesep());
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        if opt.err
            trl.cs1d3=find(trials(:,5)==4 & trials(:,8)==3 & trials(:,10)==0);
            trl.cs2d3=find(trials(:,5)==8 & trials(:,8)==3 & trials(:,10)==0);
            trl.cs1d6=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,10)==0);
            trl.cs2d6=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,10)==0);
            if min([numel(trl.cs1d3);numel(trl.cs2d3);numel(trl.cs1d6);numel(trl.cs2d6)],[],'all')<2
                continue
            end
        else
            trl.cs1d3=find(trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
            trl.cs2d3=find(trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
            trl.cs1d6=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);
            trl.cs2d6=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);
        end

        sess=['s',num2str(sessid)];

        for ff=["s1d3","s2d3","s1d6","s2d6","olf_s1","olf_s2","dur_d3","dur_d6"]
            mcid1=meta.allcid(sesssel & pct_sel.(ff));
            [~,msel1]=ismember(mcid1,suid);

            if isempty(msel1)
                continue
            end

            com_str.(['s',num2str(sessid)]).(ff).com=containers.Map('KeyType','int32','ValueType','any');
            com_str.(['s',num2str(sessid)]).(ff).com3=containers.Map('KeyType','int32','ValueType','any');
            com_str.(['s',num2str(sessid)]).(ff).com6=containers.Map('KeyType','int32','ValueType','any');
            com_str.(['s',num2str(sessid)]).(ff).com4plot=containers.Map('KeyType','int32','ValueType','any');
            
            com_str.(['s',num2str(sessid)]).(ff).fwhm=containers.Map('KeyType','int32','ValueType','any');
            com_str.(['s',num2str(sessid)]).(ff).fwhm3=containers.Map('KeyType','int32','ValueType','any');
            com_str.(['s',num2str(sessid)]).(ff).fwhm6=containers.Map('KeyType','int32','ValueType','any');

            if opt.curve % && ~opt.rnd_half
                for cc=["s1d3","s2d3","s1d6","s2d6"]
                    com_str.(['s',num2str(sessid)]).(ff).(cc)=containers.Map('KeyType','int32','ValueType','any');
                    if opt.rnd_half
                        com_str_h2.(['s',num2str(sessid)]).(ff).(cc)=containers.Map('KeyType','int32','ValueType','any');
                    end
                end
            end

            %             assert(~opt.rnd_half,'Random half-half not implemented yet')
            if opt.rnd_half

                com_str_h2.(['s',num2str(sessid)]).(ff).com=containers.Map('KeyType','int32','ValueType','any');
                com_str_h2.(['s',num2str(sessid)]).(ff).com3=containers.Map('KeyType','int32','ValueType','any');
                com_str_h2.(['s',num2str(sessid)]).(ff).com6=containers.Map('KeyType','int32','ValueType','any');
                com_str_h2.(['s',num2str(sessid)]).(ff).com4plot=containers.Map('KeyType','int32','ValueType','any');
    
                fns=fieldnames(trl);
                h1trl=struct();
                h2trl=struct();
                for fn=reshape(fns,1,[])
                    h1trl.(fn{1})=randsample(trl.(fn{1}),round(numel(trl.(fn{1}))./2));
                    h2trl.(fn{1})=trl.(fn{1})(~ismember(trl.(fn{1}),h1trl.(fn{1})));
                end
                if ismember(ff,["s1d3","s2d3","s1d6","s2d6"])
                    com_str=per_su_process(sess,suid,msel1,fr,h1trl,com_str,ff,opt);
                    com_str_h2=per_su_process(sess,suid,msel1,fr,h2trl,com_str_h2,ff,opt);
                elseif ismember(ff,["olf_s1","olf_s2"])
                    opt.currhalf=1;
                    com_str=per_su_process_olf(sess,suid,msel1,fr,h1trl,com_str,ff,opt);
                    opt.currhalf=2;
                    com_str_h2=per_su_process_olf(sess,suid,msel1,fr,h2trl,com_str_h2,ff,opt);
                else
                    com_str=per_su_process_dur(sess,suid,msel1,fr,h1trl,com_str,ff,opt);
                    com_str_h2=per_su_process_dur(sess,suid,msel1,fr,h2trl,com_str_h2,ff,opt);
                end
            else
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
            if ~opt.rnd_half && ~opt.err
                keyboard();
            end
            continue
        else
            com=sum((1:12).*mm_pref)./sum(mm_pref);
            com_str.(sess).(type).com(suid(su))=com;
            com_str.(sess).(type).fwhm(suid(su))=mm_pref2fwhm(mm_pref,12);
        end

            basemm=mean([classmm(:,1:12);classmm(3:4,13:24)],'all');
             % scale in function, gaussian smooth on plot
                SR=max(([classmm(:,1:12);classmm(3:4,13:24)]-basemm),[],'all');% removed abs
                classnnR=(classmm-basemm)./SR;
             % smooth before scale, skip smooth when plot
                if opt.early_smooth
                    class_sm=(smoothdata(classmm,2,'sgolay',5)-basemm);
                    S=max(([class_sm(:,1:12);class_sm(3:4,13:24)]),[],'all');% removed abs
                    classnn=class_sm./S;
                else
                    classnn=classnnR;
                end
                if opt.curve % && ~opt.rnd_half
                    for typeidx=1:4
                        cc=subsref(["s1d3","s2d3","s1d6","s2d6"],struct(type='()',subs={{typeidx}}));
                        if typeidx<3
                            com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,1:12);
                        else
                            com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,:);
                        end
                    end
                    % TODO: 2nd Half
                end
            if contains(type,'d3')
                com_str.(sess).(type).com4plot(suid(su))=com;
                if contains(type,'s1')
                    ppref=classnnR(1,1:12);
                else
                    ppref=classnnR(2,1:12);
                end
                ppref(ppref<0)=0;
                if any(ppref>0)
                    com=sum((1:12).*ppref)./sum(ppref);
                    com_str.(sess).(type).com3(suid(su))=com;
                    com_str.(sess).(type).fwhm3(suid(su))=mm_pref2fwhm(mm_pref,12);
                    com_str.(sess).(type).com4plot(suid(su))=com;
                end
            else
                if contains(type,'s1')
                    ppref=classnnR(3,:);
                else
                    ppref=classnnR(4,:);
                end
                ppref(ppref<0)=0;
                if any(ppref>0)
                    com=sum((1:24).*ppref)./sum(ppref);
                    com_str.(sess).(type).com6(suid(su))=com;
                    com_str.(sess).(type).fwhm6(suid(su))=mm_pref2fwhm(mm_pref,24);
                    com_str.(sess).(type).com4plot(suid(su))=com;
                end
            end
%         end
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
            if ~opt.rnd_half && ~opt.err
                keyboard();
            end
            continue
        else
            com=sum((1:12).*mm_pref)./sum(mm_pref);
            com_str.(sess).(type).com(suid(su))=com;
            com_str.(sess).(type).fwhm(suid(su))=mm_pref2fwhm(mm_pref,12);
        end
%         disp({sess,su,opt.currhalf,com})

%         if opt.curve
            basemm=mean([classmm(:,1:12);classmm(3:4,13:24)],'all');

            SR=max(([classmm(:,1:12);classmm(3:4,13:24)]-basemm),[],'all');% removed abs
            classnnR=(classmm-basemm)./SR;
         % smooth before scale, skip smooth when plot
            if opt.early_smooth
                class_sm=(smoothdata(classmm,2,'sgolay',5)-basemm);
                S=max(([class_sm(:,1:12);class_sm(3:4,13:24)]),[],'all');% removed abs
                classnn=class_sm./S;
            else
                classnn=classnnR;
            end
            if opt.curve % && ~opt.rnd_half
                for typeidx=1:4
                    cc=subsref(["s1d3","s2d3","s1d6","s2d6"],struct(type='()',subs={{typeidx}}));
                    if typeidx<3
                        com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,1:12);
                    else
                        com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,:);
                    end
                end
                % TODO: 2nd Half
            end
            if contains(type,'s1')
                ppref=classnnR(1,1:12);
            else
                ppref=classnnR(2,1:12);
            end
            ppref(ppref<0)=0;
            if any(ppref>0)
                com=sum((1:12).*ppref)./sum(ppref);
                com_str.(sess).(type).com3(suid(su))=com;
                com_str.(sess).(type).fwhm3(suid(su))=mm_pref2fwhm(mm_pref,12);
            end


            if contains(type,'s1')
                ppref=classnnR(3,:);
            else
                ppref=classnnR(4,:);
            end
            ppref(ppref<0)=0;
            if any(ppref>0)
                com=sum((1:24).*ppref)./sum(ppref);
                com_str.(sess).(type).com6(suid(su))=com;
                com_str.(sess).(type).fwhm6(suid(su))=mm_pref2fwhm(mm_pref,24);
                com_str.(sess).(type).com4plot(suid(su))=com;
            end
%         end
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
            if ~opt.rnd_half && ~opt.err
                keyboard();
            end
            continue
        else
            com=sum((1:12).*mm_pref)./sum(mm_pref);
            com_str.(sess).(type).com(suid(su))=com;
            com_str.(sess).(type).fwhm(suid(su))=mm_pref2fwhm(mm_pref,12);
        end
        
%         if opt.curve
            basemm=mean([classmm(:,1:12);classmm(3:4,13:24)],'all');
            SR=max(([classmm(:,1:12);classmm(3:4,13:24)]-basemm),[],'all');% removed abs
            classnnR=(classmm-basemm)./SR;
            if opt.early_smooth
                class_sm=(smoothdata(classmm,2,'sgolay',5)-basemm);
                S=max(([class_sm(:,1:12);class_sm(3:4,13:24)]),[],'all');% removed abs
                classnn=class_sm./S;
            else
                classnn=classnnR;
            end
            if opt.curve % && ~opt.rnd_half
                for typeidx=1:4
                    cc=subsref(["s1d3","s2d3","s1d6","s2d6"],struct(type='()',subs={{typeidx}}));
                    if typeidx<3
                        com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,1:12);
                    else
                        com_str.(sess).(type).(cc)(suid(su))=classnn(typeidx,:);
                    end
                end
                % TODO: 2nd Half
            end

            if contains(type,'d3')
                ppref=mean(classnnR(1:2,1:12));
                ppref(ppref<0)=0;
                if any(ppref>0)
                    com=sum((1:12).*ppref)./sum(ppref);
                    com_str.(sess).(type).com4plot(suid(su))=com;
                    com_str.(sess).(type).com3(suid(su))=com;
                    com_str.(sess).(type).fwhm3(suid(su))=mm_pref2fwhm(mm_pref,12);
                end
            else
                ppref=mean(classnnR(3:4,:));
                ppref(ppref<0)=0;
                if any(ppref>0)
                    com=sum((1:24).*ppref)./sum(ppref);
                    com_str.(sess).(type).com4plot(suid(su))=com;
                    com_str.(sess).(type).com6(suid(su))=com;
                    com_str.(sess).(type).fwhm6(suid(su))=mm_pref2fwhm(mm_pref,24);
                end
            end

%         end
    end
end
function fwhm=mm_pref2fwhm(mm_pref,rend)
mm_pref=mm_pref./max(mm_pref,[],'all'); % TODO: interp1
[~,xi]=max(mm_pref);
lcross=find(mm_pref(1:xi)>0.5,1);
rcross=find(mm_pref(xi:end)<0.5,1)+xi-1;
if isempty(rcross)
    rcross=rend+1;
end
fwhm=(rcross-lcross)*0.25;

end