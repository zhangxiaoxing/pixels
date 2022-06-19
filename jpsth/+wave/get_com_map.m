function com_str_=get_com_map(sel_meta,opt)
arguments
    sel_meta  (1,1) struct
    opt.onepath (1,:) char = '' % process one session under the given non-empty path
    opt.curve (1,1) logical = false % Norm. FR curve
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.rnd_half (1,1) logical = false % for bootstrap variance test
    opt.wave (1,:) char {mustBeMember(opt.wave,{'bothContext','onlyContext1','onlyContext2','anyContext1','anyContext2','any'})}  % su of wave
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} % DPA delay duration
    opt.sens_cue (1,1) double {mustBeMember(opt.sens_cue,[4,8])} % olfactory cue for duration selection
    opt.plot_COM_scheme (1,1) logical = false % for TCOM illustration
    opt.one_SU_showcase (1,1) logical = false % for the TCOM-FC joint showcase
    opt.append_late_delay (1,1) logical = false % Uses stats from early delay but include illustration for late delay

end
persistent com_str opt_ sel_meta_

if isempty(com_str) || ~isequaln(opt,opt_) || ~isequaln(sel_meta_,sel_meta) 
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

    %% RANKSUM
    waveid=sel_meta.wave_id;
    % 3onlyS1:1,3onlyS2:2,6onlyS1:3,6onlyS2:4,bothS1:5,bothS2:6,nonmem:0,switch,etc:-1,
    switch opt.wave
        case 'bothContext' %both 3 and 6s
            c1ids=5;c2ids=6;
        case 'onlyContext1'
            c1ids=1;c2ids=2;
        case 'onlyContext2'
            c1ids=3;c2ids=4;
        case 'anyContext1'
            c1ids=[1 5];c2ids=[2 6];
        case 'anyContext2'
            c1ids=[3 5];c2ids=[4 6];
        case 'any'
            c1ids=[1 3 5];c2ids=[2 4 6];
    end

    %%
    for sessid=reshape(usess,1,[])
        fpath=fullfile(homedir,ephys.sessid2path(sessid),'FR_All_ 250.hdf5');
        sesssel=meta.sess==sessid;
        if ~any(sesssel), continue;end
        fpath=replace(fpath,'\',filesep());
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        %% RANKSUM
%         mcid1=meta.allcid(ismember(waveid,c1ids) & sesssel & ismember(meta.reg_tree(1,:),{'CH','BS'}).');
%         mcid2=meta.allcid(ismember(waveid,c2ids) & sesssel & ismember(meta.reg_tree(1,:),{'CH','BS'}).');
        mcid1=meta.allcid(ismember(waveid,c1ids) & sesssel);
        mcid2=meta.allcid(ismember(waveid,c2ids) & sesssel);

        msel1=find(ismember(suid,mcid1));
        msel2=find(ismember(suid,mcid2));

        if opt.rnd_half
            for ff=["c1a","c2a","c1b","c2b","c1e","c2e"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["c1aheat","c2aheat","c1acurve","c2acurve","c1aanticurve","c2aanticurve",...
                        "c1bheat","c2bheat","c1bcurve","c2bcurve","c1banticurve","c2banticurve",...
                        "c1eheat","c2eheat","c1ecurve","c2ecurve","c1eanticurve","c2eanticurve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
        else
            for ff=["c1","c2"]
                com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
            end
            if opt.curve
                for ff=["c1heat","c2heat","c1curve","c2curve","c1anticurve","c2anticurve"]
                    com_str.(['s',num2str(sessid)]).(ff)=containers.Map('KeyType','int32','ValueType','any');
                end
            end
        end

        if isempty(msel1) && isempty(msel2)
            continue
        end
        
        %% TODO Needs to handle duration groups
        if isfield(sel_meta,'selec_d3') % Sensory wave 
            opt.dur=false;
            c1sel=find(trials(:,5)==4 & trials(:,8)==opt.delay & trials(:,9)>0 & trials(:,10)>0);
            c2sel=find(trials(:,5)==8 & trials(:,8)==opt.delay & trials(:,9)>0 & trials(:,10)>0);

            e1sel=find(trials(:,5)==4 & trials(:,8)==opt.delay & trials(:,10)==0);
            e2sel=find(trials(:,5)==8 & trials(:,8)==opt.delay & trials(:,10)==0);
        elseif isfield(sel_meta,'selec_s1') % Duration wave
            opt.dur=true;
            opt.append_late_delay=true;
            opt.delay=3; % this does truncate the wave heatmap plot
            c1sel=find(trials(:,5)==opt.sens_cue & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
            c2sel=find(trials(:,5)==opt.sens_cue & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);

            e1sel=find(trials(:,5)==opt.sens_cue & trials(:,8)==3 & trials(:,10)==0);
            e2sel=find(trials(:,5)==opt.sens_cue & trials(:,8)==6 & trials(:,10)==0);
        else
            error('Data mismatch')
        end
        sess=['s',num2str(sessid)];

        if opt.rnd_half
            c1a=randsample(c1sel,floor(numel(c1sel)./2));
            %             c1a=c1sel(1:2:end);
            c1b=c1sel(~ismember(c1sel,c1a));

            c2a=randsample(c2sel,floor(numel(c2sel)./2));
            %             c2a=c2sel(1:2:end);
            c2b=c2sel(~ismember(c2sel,c2a));
            if nnz(c1a)>2 && nnz(c1b)>2 && nnz(c2a)>2 && nnz(c2b)>2 && nnz(e1sel)>2 && nnz(e2sel)>2
                com_str=per_su_process(sess,suid,msel1,fr,c1a,c2a,com_str,'c1a',opt);
                com_str=per_su_process(sess,suid,msel1,fr,c1b,c2b,com_str,'c1b',opt);
                com_str=per_su_process(sess,suid,msel2,fr,c2a,c1a,com_str,'c2a',opt);
                com_str=per_su_process(sess,suid,msel2,fr,c2b,c1b,com_str,'c2b',opt);
                com_str=per_su_process(sess,suid,msel1,fr,e1sel,e2sel,com_str,'c1e',opt);
                com_str=per_su_process(sess,suid,msel2,fr,e2sel,e1sel,com_str,'c2e',opt);
            end
        else
            com_str=per_su_process(sess,suid,msel1,fr,c1sel,c2sel,com_str,'c1',opt);
            com_str=per_su_process(sess,suid,msel2,fr,c2sel,c1sel,com_str,'c2',opt);
        end
    end
end
com_str_=com_str;
opt_=opt;
sel_meta_=sel_meta;
end

function com_str=per_su_process(sess,suid,msel,fr,pref_sel,nonpref_sel,com_str,samp,opt)
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

        mm_pref=squeeze(mean(fr(pref_sel,su,stats_window))).'-basemm;

        if opt.append_late_delay
            if strcmp(samp,'c1') %i.e. 3s-duration
                curve=mm_pref;
                anticurve=squeeze(mean(fr(nonpref_sel,su,17:40))).'-basemm;
            else % 'c2', i.e. 6s-duration
                curve=squeeze(mean(fr(pref_sel,su,17:40))).'-basemm;
                anticurve=squeeze(mean(fr(nonpref_sel,su,stats_window))).'-basemm;
            end
        else
            curve=mm_pref;
            anticurve=squeeze(mean(fr(nonpref_sel,su,stats_window))).'-basemm;
        end
        mm_pref(mm_pref<0)=0;

        if ~any(mm_pref>0)
            if (ismember(samp,{'c1e','c2e'})) ...
                    || (strcmp(opt.wave,'onlyContext2') && opt.delay==3 && evalin('caller','isfield(sel_meta,''selec_d3'')'))... % TODO: fix recursive call stack
                    || (strcmp(opt.wave,'onlyContext1') && opt.delay==6 && evalin('caller','isfield(sel_meta,''selec_d3'')')) % TODO: fix recursive call stack
                com=-1;
                disp(string(samp)+" PEAK mismatch, TCOM set to -1")
            else
                disp(string(samp)+" NOPEAK")
                keyboard();
            end
        else
            com=sum((1:numel(stats_window)).*mm_pref)./sum(mm_pref);
        end
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
                    curve=curve(1:numel(stats_window));
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

