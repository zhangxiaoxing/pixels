function com_str_=get_com_map(opt)
arguments
    opt.sel_meta  (1,1) struct
    opt.onepath (1,:) char = '' % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
    opt.cell_type (1,:) char {mustBeMember(opt.cell_type,{'any_s1','any_s2','grey_sel'})} = 'grey_sel' % select (sub) populations
    opt.selidx (1,1) logical = false % calculate COM of selectivity index
    opt.wave (1,:) char {mustBeMember(opt.wave,{'both','only3','only6','any3','any6','any'})}  % su of wave
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} % DPA delay duration
    opt.plot_COM_scheme (1,1) logical = false % for TCOM illustration
    opt.one_SU_showcase (1,1) logical = false % for the TCOM-FC joint showcase
    opt.stats_model (1,:) char {mustBeMember(opt.stats_model,{'ANOVA2','ANOVA3','RANKSUM','RANKSUM2'})} = 'RANKSUM'  % 2-way anova, 3-way anova (with time bin), s1-s2 rannksum, 2-way ranksum (w/ duration)

end
persistent com_str opt_
wrs_stats=strcmp(opt.stats_model,'RANKSUM') & isfield(opt,'sel_meta') & isfield(opt.sel_meta,'wave_id');
assert(wrs_stats);

if true || isempty(com_str) || ~isequaln(opt,opt_)
    meta=ephys.util.load_meta();

    homedir=ephys.util.getHomedir('type','raw');
    fl=dir(fullfile(homedir,'**','FR_All_ 250.hdf5'));
    com_str=struct();

    %% RANKSUM
    if strcmp(opt.stats_model,'RANKSUM')
        waveid=opt.sel_meta.wave_id;
        % 3onlyS1:1,3onlyS2:2,6onlyS1:3,6onlyS2:4,bothS1:5,bothS2:6,nonmem:0,switch,etc:-1,
        switch opt.wave
            case 'both' %both 3 and 6s
                s1ids=5;s2ids=6;
            case 'only3'
                s1ids=1;s2ids=2;
            case 'only6'
                s1ids=3;s2ids=4;
            case 'any3'
                s1ids=[1 5];s2ids=[2 6];
            case 'any6'
                s1ids=[3 5];s2ids=[4 6];
            case 'any'
                s1ids=[1 3 5];s2ids=[2 4 6];
        end
%TODO: ANOVA2        
%     elseif strcmp(opt.stats_model,'ANOVA2')
%         switch opt.wave %both 3 and 6s
%             case 'both'
%                 s1ids=5;s2ids=6;
%             case 'only3'
%                 s1ids=1;s2ids=2;
%             case 'only6'
%                 s1ids=3;s2ids=4;
%             case 'any3'
%                 s1ids=[1 5];s2ids=[2 6];
%             case 'any6'
%                 s1ids=[3 5];s2ids=[4 6];
%             case 'any'
%                 s1ids=[1 3 5];s2ids=[2 4 6];
%         end
    end

    %%
    for ii=1:size(fl,1)
        if strlength(opt.onepath)==0
            dpath=regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once');
            fpath=fullfile(fl(ii).folder,fl(ii).name);
        else
            dpath=regexp(opt.onepath,'(?<=SPKINFO[\\/]).*$','match','once');
            if isempty(dpath)
                dpath=opt.onepath;
            end
            fpath=fullfile(homedir,dpath,'FR_All_ 250.hdf5');
        end
        fpath=replace(fpath,'\',filesep());
        pc_stem=replace(dpath,'/','\');
        sesssel=startsWith(meta.allpath,pc_stem);
        if ~any(sesssel), continue;end
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        %% RANKSUM
        if strcmp(opt.stats_model,'RANKSUM')
            switch opt.cell_type
                case 'any_s1'
                    mcid1=meta.allcid(ismember(waveid,s1ids) & sesssel);
                    mcid2=[];
                case 'any_s2'
                    mcid1=[];
                    mcid2=meta.allcid(ismember(waveid,s2ids) & sesssel);
                case 'grey_sel'
                    mcid1=meta.allcid(ismember(waveid,s1ids) & sesssel & ismember(meta.reg_tree(1,:),{'CH','BS'}).');
                    mcid2=meta.allcid(ismember(waveid,s2ids) & sesssel & ismember(meta.reg_tree(1,:),{'CH','BS'}).');
            end
        elseif strcmp(opt.stats_model,'ANOVA2')
            
        end
        %%

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
        %% TODO Needs to handle duration groups
        s1sel=find(trials(:,5)==4 & trials(:,8)==opt.delay & trials(:,9)>0 & trials(:,10)>0);
        s2sel=find(trials(:,5)==8 & trials(:,8)==opt.delay & trials(:,9)>0 & trials(:,10)>0);

        e1sel=find(trials(:,5)==4 & trials(:,8)==opt.delay & trials(:,10)==0);
        e2sel=find(trials(:,5)==8 & trials(:,8)==opt.delay & trials(:,10)==0);

        sess=['s',num2str(sessid)];
        %     if sum(trial(:,9))<40,continue;end %meta data obtained from processed
        %     welltrained dataset
        if opt.rnd_half
            for ff=["s1a","s2a","s1b","s2b","s1e","s2e"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["s1aheat","s2aheat","s1acurve","s2acurve","s1aanticurve","s2aanticurve",...
                        "s1bheat","s2bheat","s1bcurve","s2bcurve","s1banticurve","s2banticurve",...
                        "s1eheat","s2eheat","s1ecurve","s2ecurve","s1eanticurve","s2eanticurve"]
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
                for ff=["s1heat","s2heat","s1curve","s2curve","s1anticurve","s2anticurve"]
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
opt_=opt;
end

function com_str=per_su_process(sess,suid,msel,fr,pref_sel,nonpref_sel,com_str,samp,opt)
%TODO if decision
if opt.decision
    stats_window=(opt.delay*4+17):44;
else
    stats_window=17:(opt.delay*4+16);
end
for su=reshape(msel,1,[])
    perfmat=squeeze(fr(pref_sel,su,:));
    npmat=squeeze(fr(nonpref_sel,su,:));
    basemm=mean([mean(perfmat(:,stats_window),1);mean(npmat(:,stats_window),1)]);

    basemm=mean(basemm);

    itimm=mean(fr([pref_sel;nonpref_sel],su,1:12),'all');
    %TODO compare the effect of smooth
    if strcmp(opt.cell_type,'any_nonmem')
        mm=smooth(squeeze(mean(fr([pref_sel;nonpref_sel],su,:))),5).';
        mm_pref=mm(stats_window)-itimm;
    else
        mm=smooth(squeeze(mean(fr(pref_sel,su,:))),5).';
        mm_pref=mm(stats_window)-basemm;
    end
    if max(mm_pref)<=0,disp('NOPEAK');continue;end % work around 6s paritial
    if opt.selidx
        sel_vec=[mean(perfmat,1);mean(npmat,1)];
        sel_idx=(-diff(sel_vec)./sum(sel_vec));
        sel_idx(all(sel_vec==0))=0;
        curve=sel_idx(stats_window);
    elseif contains(opt.cell_type,'any')
        curve=squeeze(mean(fr(pref_sel,su,:))).'-itimm;
        anticurve=squeeze(mean(fr(nonpref_sel,su,:))).'-itimm;
    elseif opt.delay==6 % full, early or late
        curve=squeeze(mean(fr(pref_sel,su,17:40))).'-basemm;
        anticurve=squeeze(mean(fr(nonpref_sel,su,17:40))).'-basemm;
    else
        curve=mm_pref;
        anticurve=squeeze(mean(fr(nonpref_sel,su,stats_window))).'-basemm;
    end
    mm_pref(mm_pref<0)=0;
    com=sum((1:numel(stats_window)).*mm_pref)./sum(mm_pref);
    if opt.plot_COM_scheme
        fc_scheme(curve,mm_pref,com)
    end
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

%% for COM scheme illustration
function fc_scheme(curve,mm_pref,com)
curve(curve<0)=0;
close all
fh=figure('Color','w', 'Position',[32,32,275,225]);
bar(mm_pref,'k');
ylabel('Baseline-subtracted firing rate w (Hz)')
set(gca,'XTick',0:4:24,'XTickLabel',0:6)
xlabel('Time t (0.25 to 6 sec in step of 0.25 sec)')
xline(com,'--r');
keyboard()
%             exportgraphics(gcf(),'COM_illustration_L0.pdf','ContentType','vector')

end

