%TODO Error trial
function com_str_=get_com_map(opt)
arguments
    opt.onepath (1,:) char = '' % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.per_sec_stats (1,1) logical = false % calculate COM using per-second mean as basis for normalized firing rate, default is coss-delay mean
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
    opt.keep_sust (1,1) logical = false % use sustained coding neuron
    opt.selidx (1,1) logical = false % calculate COM of selectivity index
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6 % DPA delay duration
end
persistent com_str onepath_ delay_ selidx_ decision_ rnd_half_ curve_

if opt.delay==6
    warning('Delay set to default 6')
end

if isempty(onepath_), onepath_='';end
if isempty(com_str) || ~strcmp(opt.onepath, onepath_) || opt.delay~=delay_ || opt.selidx~=selidx_ || opt.decision~=decision_ || opt.rnd_half~=rnd_half_ || opt.curve~=curve_
    meta_str=ephys.util.load_meta('type','neupix','delay',opt.delay);
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
        if opt.keep_sust
            mcid1=meta_str.allcid(ismember(meta_str.mem_type,1:2) & sesssel.'& strcmp(meta_str.reg_tree(2,:),'CTX'));
            mcid2=meta_str.allcid(ismember(meta_str.mem_type,3:4) & sesssel.'& strcmp(meta_str.reg_tree(2,:),'CTX'));
        else
            mcid1=meta_str.allcid(meta_str.mem_type==2 & sesssel.' & strcmp(meta_str.reg_tree(2,:),'CTX'));
            mcid2=meta_str.allcid(meta_str.mem_type==4 & sesssel.' & strcmp(meta_str.reg_tree(2,:),'CTX'));
        end
        msel1=find(ismember(suid,mcid1));
        msel2=find(ismember(suid,mcid2));
        if isempty(msel1) && isempty(msel2)
            if strlength(opt.onepath)==0
                continue
            else
                break
            end
        end
        sessid=ephys.path2sessid(pc_stem);
        s1sel=find(trial(:,5)==4 & trial(:,8)==opt.delay & trial(:,9)>0 & trial(:,10)>0);
        s2sel=find(trial(:,5)==8 & trial(:,8)==opt.delay & trial(:,9)>0 & trial(:,10)>0);
        
        e1sel=find(trial(:,5)==4 & trial(:,8)==opt.delay & trial(:,10)==0);
        e2sel=find(trial(:,5)==8 & trial(:,8)==opt.delay & trial(:,10)==0);
        
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
                        "s1eheat","s2eheat","s1ecurve","s2ecurve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
            s1a=randsample(s1sel,floor(numel(s1sel)./2));
%             s1a=s1sel(1:2:end);
            s1b=s1sel(~ismember(s1sel,s1a));

            s2a=randsample(s2sel,floor(numel(s2sel)./2));
%             s2a=s2sel(1:2:end);
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
delay_ =opt.delay;
selidx_=opt.selidx;
decision_=opt.decision;
rnd_half_=opt.rnd_half;
curve_=opt.curve;
end

function com_str=per_su_process(sess,suid,msel,fr,pref_sel,nonpref_sel,com_str,samp,opt)
%TODO if decision
if opt.decision
    stats_window=(opt.delay*4+17):44;
else
    stats_window=17:(opt.delay*4+16);
end
for su=reshape(msel,1,[])
    perfmat=squeeze(fr(pref_sel,su,stats_window));
    npmat=squeeze(fr(nonpref_sel,su,stats_window));
    basemm=mean([mean(perfmat,1);mean(npmat,1)]);
    if opt.selidx
        sel_vec=([mean(perfmat,1);mean(npmat,1)]);
        sel_idx=(-diff(sel_vec)./sum(sel_vec));
        curve=sel_idx;
        sel_idx(sel_idx<0 | all(sel_vec==0))=0;
        com=sum((1:numel(stats_window)).*sel_idx)./sum(sel_idx);
        if ~isfinite(com)
            fprintf('Moved sess %d su%d, infinite TCOM\n',sess,su)
            com=numel(stats_window)+1;
        end
    else
        if ~opt.per_sec_stats
            basemm=mean(basemm);
        end
        %TODO check the effect of smooth
        mm=smooth(squeeze(mean(fr(pref_sel,su,:))),5).';
        mm_pref=mm(stats_window)-basemm;
        if max(mm_pref)<=0,continue;end
        curve=mm_pref;
        mm_pref(mm_pref<0)=0;
        com=sum((1:numel(stats_window)).*mm_pref)./sum(mm_pref);

        %% for COM showcase
%         if min(curve)>0
% %             close all
%             fh=figure('Color','w');
%             bar(mm_pref,'k');
%             ylabel('Baseline-deduced firing rate w (Hz)')
%             set(gca,'XTick',0:4:24,'XTickLabel',0:6)
%             xlabel('Time t (0.25 to 6 sec in step of 0.25 sec)')
%             xline(com,'--r');
%             keyboard()
% %             exportgraphics(gcf(),'COM_showcase.pdf','ContentType','vector')
%         end

    end
    com_str.(sess).(samp)(suid(su))=com;
    if opt.curve
        com_str.(sess).([samp,'curve'])(suid(su))=curve;
        if opt.rnd_half
            heatnorm=curve./max(curve);
        else % per_su_showcase
            heatcent=squeeze(fr(pref_sel,su,stats_window))-basemm; %centralized norm. firing rate for heatmap plot
            heatnorm=heatcent./max(abs(heatcent));
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
