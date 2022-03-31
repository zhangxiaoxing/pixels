function com_str_=get_dur_com_map(opt)
arguments
    opt.onepath (1,:) char = '' % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
end
persistent com_str opt_

if true || isempty(com_str) || ~isequaln(opt,opt_)
    meta=ephys.util.load_meta();
    anovameta=wave.get_dur_waveid();

    homedir=ephys.util.getHomedir('type','raw');

    com_str=struct();
    usess=unique(anovameta.sess);
    for sessid=reshape(usess,1,[])
        sesssel=anovameta.sess==sessid;
        if ~any(sesssel), continue;end

        if strlength(opt.onepath)==0
            fpath=fullfile(homedir,ephys.sessid2path(sessid),'FR_All_ 250.hdf5');
        else
            dpath=regexp(opt.onepath,'(?<=SPKINFO[\\/]).*$','match','once');
            if isempty(dpath)
                dpath=opt.onepath;
            end
            fpath=fullfile(homedir,dpath,'FR_All_ 250.hdf5');
        end
        fpath=replace(fpath,'\',filesep());
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        mcid1=anovameta.allcid(anovameta.dur_waveid==3 & sesssel & ismember(meta.reg_tree(1,:),{'CH','BS'}).');
        mcid2=anovameta.allcid(anovameta.dur_waveid==6 & sesssel & ismember(meta.reg_tree(1,:),{'CH','BS'}).');

        msel1=find(ismember(suid,mcid1));
        msel2=find(ismember(suid,mcid2));
        if isempty(msel1) && isempty(msel2) && strlength(opt.onepath)~=0
            break
        end

        d3sel=find(ismember(trials(:,5),[4,8]) & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
        d6sel=find(ismember(trials(:,5),[4 8]) & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);

        e3sel=find(ismember(trials(:,5),[4,8]) & trials(:,8)==3 & trials(:,10)==0);
        e6sel=find(ismember(trials(:,5),[4 8]) & trials(:,8)==6 & trials(:,10)==0);
        sess=['s',num2str(sessid)];
        if opt.rnd_half
            for ff=["d3a","d6a","d3b","d6b","d3e","d6e"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["d3aheat","d6aheat","d3acurve","d6acurve","d3aanticurve","d6aanticurve",...
                        "d3bheat","d6bheat","d3bcurve","d6bcurve","d3banticurve","d6banticurve",...
                        "d3eheat","d6eheat","d3ecurve","d6ecurve","d3eanticurve","d6eanticurve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
            d3a=randsample(d3sel,floor(numel(d3sel)./2));
            d3b=d3sel(~ismember(d3sel,d3a));

            d6a=randsample(d6sel,floor(numel(d6sel)./2));
            d6b=d6sel(~ismember(d6sel,d6a));
            if nnz(d3a)>2 && nnz(d3b)>2 && nnz(d6a)>2 && nnz(d6b)>2 && nnz(e3sel)>2 && nnz(e6sel)>2
                com_str=per_su_process(sess,suid,msel1,fr,d3a,d6a,com_str,'d3a',opt);
                com_str=per_su_process(sess,suid,msel1,fr,d3b,d6b,com_str,'d3b',opt);
                com_str=per_su_process(sess,suid,msel2,fr,d6a,d3a,com_str,'d6a',opt);
                com_str=per_su_process(sess,suid,msel2,fr,d6b,d3b,com_str,'d6b',opt);
                com_str=per_su_process(sess,suid,msel1,fr,e3sel,e6sel,com_str,'d3e',opt);
                com_str=per_su_process(sess,suid,msel2,fr,e6sel,e3sel,com_str,'d6e',opt);
            end
        else
            for ff=["d3","d6"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["d3heat","d6heat","d3curve","d6curve","d3anticurve","d6anticurve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
            com_str=per_su_process(sess,suid,msel1,fr,d3sel,d6sel,com_str,'d3',opt);
            com_str=per_su_process(sess,suid,msel2,fr,d6sel,d3sel,com_str,'d6',opt);
        end
        if ~strlength(opt.onepath)==0
            break;
        end
    end
end
com_str_=com_str;
opt_=opt;
end

function com_str=per_su_process(sess,suid,msel,fr,pref_sel,nonpref_sel,com_str,samp,opt)
stats_window=17:28;
for su=reshape(msel,1,[])
    perfmat=squeeze(fr(pref_sel,su,:));
    npmat=squeeze(fr(nonpref_sel,su,:));
    basemm=mean([mean(perfmat(:,stats_window),1);mean(npmat(:,stats_window),1)]);

    basemm=mean(basemm);
    
    %TODO compare the effect of smooth
    mm=smooth(squeeze(mean(fr(pref_sel,su,:))),5).';
    mm_pref=mm(stats_window)-basemm;

    if max(mm_pref)<=0,disp('NOPEAK');continue;end % work around 6s paritial
    curve=mm_pref;
    anticurve=squeeze(mean(fr(nonpref_sel,su,stats_window))).'-basemm;
    mm_pref(mm_pref<0)=0;
    com=sum((1:numel(stats_window)).*mm_pref)./sum(mm_pref);
    com_str.(sess).(samp)(suid(su))=com;
    if opt.curve
        com_str.(sess).([samp,'curve'])(suid(su))=curve;
        if exist('anticurve','var')
            com_str.(sess).([samp,'anticurve'])(suid(su))=anticurve;
        end
        heatcent=squeeze(fr(pref_sel,su,stats_window))-basemm; %centralized norm. firing rate for heatmap plot
        heatnorm=heatcent./max(abs(heatcent));
        heatnorm(heatnorm<0)=0;
        if size(heatnorm,1)>10
            if numel(curve)>numel(stats_window)
                curve=curve(stats_window);
            end
            cc=arrayfun(@(x) min(corrcoef(heatnorm(x,:),curve),[],'all'),1:size(heatnorm,1));
            [~,idx]=sort(cc,'descend');
            heatnorm=heatnorm(idx(1:10),:);
        end
        com_str.(sess).([samp,'heat'])(suid(su))=heatnorm;
    end
end
end

