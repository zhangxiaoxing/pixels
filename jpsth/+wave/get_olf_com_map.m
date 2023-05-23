function [com_str_,com_str_h2]=get_olf_com_map(sel_meta,opt)
arguments
    sel_meta  (1,1) struct
    opt.one_sess (1,:) double = [] % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
    opt.err (1,1) logical = false %u stats in error trials
end

%TODO proper declaration
persistent com_str opt_ sel_meta_

assert(~(opt.rnd_half && opt.err), 'Unable to handle 2Fold CV and error trial in same run')
if opt.rnd_half || isempty(com_str) || ~isequaln(opt,opt_) || ~isequaln(sel_meta_,sel_meta) 
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

    wm_sel=struct();

    wm_sel.olf_s1=sel_meta.wave_id==5;
    wm_sel.olf_s2=sel_meta.wave_id==6;

    for sessid=reshape(usess,1,[])
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

        for ff=["olf_s1","olf_s2"]
            mcid1=meta.allcid(sesssel & wm_sel.(ff));
            [~,msel1]=ismember(mcid1,suid);
            if isempty(msel1)
                continue
            end

            com_str.(['s',num2str(sessid)]).(ff).com=containers.Map('KeyType','int32','ValueType','any');
            com_str.(['s',num2str(sessid)]).(ff).fwhm=containers.Map('KeyType','int32','ValueType','any');

            if opt.curve
                for cc=["s1","s2"]
                    com_str.(['s',num2str(sessid)]).(ff).(cc)=containers.Map('KeyType','int32','ValueType','any');
                    if opt.rnd_half
                        com_str_h2.(['s',num2str(sessid)]).(ff).(cc)=containers.Map('KeyType','int32','ValueType','any');
                    end
                end
            end

            if opt.rnd_half
                com_str_h2.(['s',num2str(sessid)]).(ff).com=containers.Map('KeyType','int32','ValueType','any');
                fns=fieldnames(trl);
                h1trl=struct();
                h2trl=struct();
                for fn=reshape(fns,1,[])
                    h1trl.(fn{1})=randsample(trl.(fn{1}),round(numel(trl.(fn{1}))./2));
                    h2trl.(fn{1})=trl.(fn{1})(~ismember(trl.(fn{1}),h1trl.(fn{1})));
                end
                com_str=per_su_process_olf(sess,suid,msel1,fr,h1trl,com_str,ff,opt);
                com_str_h2=per_su_process_olf(sess,suid,msel1,fr,h2trl,com_str_h2,ff,opt);
            else
                com_str=per_su_process_olf(sess,suid,msel1,fr,trl,com_str,ff,opt);
            end
        end
    end
end
com_str_=com_str;
opt_=opt;
sel_meta_=sel_meta;
end

function com_str=per_su_process_olf(sess,suid,msel,fr,trls,com_str,type,opt)
    for su=reshape(msel,1,[])
        classmm=[mean([squeeze(fr(trls.cs1d3,su,17:28));squeeze(fr(trls.cs1d6,su,17:28))]),... %s1d3
            mean([squeeze(fr(trls.cs1d6,su,29:40))]);...%s1d6
            mean([squeeze(fr(trls.cs2d3,su,17:28));squeeze(fr(trls.cs2d6,su,17:28))]),...%s2d3
            mean([squeeze(fr(trls.cs2d6,su,29:40))])...%s2d6
            ];
        basemm=mean(classmm,'all');
        
        S=max(classmm-basemm,[],'all');% removed abs
        classnn=(classmm-basemm)./S;
        
        if contains(type,'s1')
            mm_pref=classnn(1,:);
        else
            mm_pref=classnn(2,:);
        end
        mm_pref(mm_pref<0)=0;
        if ~any(mm_pref>0)
            if ~opt.rnd_half && ~opt.err
                keyboard();
            else
                continue
            end
        else
            com=sum((1:24).*mm_pref)./sum(mm_pref);
            com_str.(sess).(type).com(suid(su))=com;
            com_str.(sess).(type).fwhm(suid(su))=mm_pref2fwhm(mm_pref,24);
        end
        % smooth early option removed due to necessity of smooth along
        % the neuron# axis, due to limited pixels in printing or display
        if opt.curve
            com_str.(sess).(type).s1(suid(su))=classnn(1,:);
            com_str.(sess).(type).s2(suid(su))=classnn(2,:);
        end
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