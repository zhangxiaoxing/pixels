if ~exist('fidx','var') || ~isfile(fullfile('.',sprintf('0203_BZ_XCORR_duo_f%d.mat',fidx)))
    disp('Error loading file')
    if isunix
        quit(0);
    else
        return
    end
end
disp(fidx)
load(fullfile('.',sprintf('0203_BZ_XCORR_duo_f%d.mat',fidx)))
% ts_sep=0.2;
ts_sep=0;
pre_thresh=0;
trl_thresh=5;
%%
fc.util.loaded
fc.util.dependency
%%
sums=cell(0);

for idx=1:size(mono.sig_con,1)
    trials=h5read(fullfile(folder,'events.hdf5'),'/trials')';
    ctrials=behav.procPerf(trials);
    fstr=load(fullfile(folder,'spike_info.mat'));
    spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
    spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);
    FT_SPIKE=struct();
    cluster_ids=mono.sig_con(idx,:);
    FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids')));
    FT_SPIKE.timestamp=cell(1,2);
    for i=1:numel(cluster_ids)
        FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
    end
    %  continuous format F T struct file
    sps=30000;
    cfg=struct();
    cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
    cfg.trlunit='timestamps';
    cfg.timestampspersecond=sps;
    FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
    
    %dimord={ts1,ts2, fcevents, anti-causal fc, 50ms-lag fc} x trl(~240) x bin(-3:11)
    
    stats=zeros(5,size(trials,1),14);
    for trlIdx=1:size(trials,1)
        ts1=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==trlIdx)';
        ts2=FT_SPIKE.time{2}(FT_SPIKE.trial{2}==trlIdx)';
        ts1=ts1([false;diff(ts1)>ts_sep]);
        if isempty(ts1)
            continue
        end
        ts1(:,2)=1;
        ts2(:,2)=2;
        ts_id=[ts1;ts2];
        [~,s]=sort(ts_id(:,1));
        ts_id=ts_id(s,:);
        ts_id_tagged=fc.fc_tag(ts_id,false);
        stats(1,trlIdx,:)=histcounts(ts1(:,1),-3:11);
        stats(2,trlIdx,:)=histcounts(ts2(:,1),-3:11);
        stats(3,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,3)==1 & ts_id_tagged(:,2)==1,1),-3:11);
        stats(4,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,4)==1 & ts_id_tagged(:,2)==1,1),-3:11);
        stats(5,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,5)==1 & ts_id_tagged(:,2)==1,1),-3:11);
    end
    
    sel_6s_S1=trials(:,5)==4 & trials(:,8)== 6;
    sel_6s_S2=trials(:,5)==8 & trials(:,8)== 6;
    
    sel_3s_S1=trials(:,5)==4 & trials(:,8)== 3;
    sel_3s_S2=trials(:,5)==8 & trials(:,8)== 3;
    
    %pair x bin x {mean|S1, mean|S2, n|S1, n|S2, wrs-p} x {3s, 6s}:
    onepair=nan(14,21,2);
    for sbin=1:14
        pre_sel=stats(1,:,sbin)'>pre_thresh;
        s31fc=stats(3,sel_3s_S1 & pre_sel,sbin);
        s31pre=stats(1,sel_3s_S1 & pre_sel,sbin);
        s31post=stats(2,sel_3s_S1 & pre_sel,sbin);

        s32fc=stats(3,sel_3s_S2 & pre_sel,sbin);
        s32pre=stats(1,sel_3s_S2 & pre_sel,sbin);
        s32post=stats(2,sel_3s_S1 & pre_sel,sbin);

        s61fc=stats(3,sel_6s_S1 & pre_sel,sbin);
        s61pre=stats(1,sel_6s_S1 & pre_sel,sbin);
        s61post=stats(2,sel_6s_S1 & pre_sel,sbin);
        
        s62fc=stats(3,sel_6s_S2 & pre_sel,sbin);
        s62pre=stats(1,sel_6s_S2 & pre_sel,sbin);
        s62post=stats(2,sel_6s_S2 & pre_sel,sbin);
        
        
        if nnz(sel_3s_S1 & pre_sel)>=trl_thresh && nnz(sel_3s_S2 & pre_sel)>=trl_thresh && ...
                nnz(sel_6s_S1 & pre_sel)>=trl_thresh && nnz(sel_6s_S2 & pre_sel)>=trl_thresh
%             s31frac=s31fc./s31pre;
%             s32frac=s32fc./s32pre;
%             s61frac=s61fc./s61pre;
%             s62frac=s62fc./s62pre;
%             if nnz([s31frac,s32frac])==0
%                 p3=1;
%             else
%                 p3=ranksum(s31frac,s32frac);
%             end
%             if nnz([s61frac,s62frac])==0
%                 p6=1;
%             else
%                 p6=ranksum(s61frac,s62frac);
%             end
            
%             onepair(sbin,:,1)=[mean(s31frac),mean(s32frac),sum(s31pre),sum(s32pre),p3];
%             onepair(sbin,:,2)=[mean(s61frac),mean(s62frac),sum(s61pre),sum(s62pre),p6];

            if nnz([s31pre,s32pre])==0
                p3pre=1;
            else
                p3pre=ranksum(s31pre,s32pre);
            end
            if nnz([s61pre,s62pre])==0
                p6pre=1;
            else
                p6pre=ranksum(s61pre,s62pre);
            end
            

            if nnz([s31post,s32post])==0
                p3post=1;
            else
                p3post=ranksum(s31post,s32post);
            end
            if nnz([s61post,s62post])==0
                p6post=1;
            else
                p6post=ranksum(s61post,s62post);
            end            
            

            if nnz([s31fc,s32fc])==0
                p3=1;
            else
                p3=ranksum(s31fc,s32fc);
            end
            if nnz([s61fc,s62fc])==0
                p6=1;
            else
                p6=ranksum(s61fc,s62fc);
            end


            onepair(sbin,:,1)=[mean(s31fc),std(s31fc),numel(s31fc),...
                mean(s32fc),std(s32fc),numel(s32fc),...
                mean(s31pre),std(s31pre),numel(s31pre),...
                mean(s32pre),std(s32pre),numel(s32pre),...
                mean(s31post),std(s31post),numel(s31post),...
                mean(s32post),std(s32post),numel(s32post),...
                p3,p3pre,p3post];

                onepair(sbin,:,2)=[mean(s61fc),std(s61fc),numel(s61fc),...
                mean(s62fc),std(s62fc),numel(s62fc),...
                mean(s61pre),std(s61pre),numel(s61pre),...
                mean(s62pre),std(s62pre),numel(s62pre),...
                mean(s61post),std(s61post),numel(s61post),...
                mean(s62post),std(s62post),numel(s62post),...
                p6,p6pre,p6post];


        else
            continue
        end
    end
    sums(end+1,:)={idx,mono.sig_con(idx,:),onepair};
    if exist('debug','var') && debug && size(sums,1)>=50
        break
    end
end
save(sprintf('fc_coding_f%d.mat',fidx),'sums')
if isunix
    quit(0)
else
    return
end
% 
% for i=1:size(sums,1)
%     if nnz([sums{i,3}(:,5,1);sums{i,3}(:,5,2)]<0.01)>2
% %         disp(sums{i,3})
%         fh=figure();
%         hold on
%         plot((-2:11)-0.5,sums{i,3}(:,1,1),'m-');
%         plot((-2:11)-0.5,sums{i,3}(:,2,1),'c-');
%         plot((-2:11)-0.5,sums{i,3}(:,1,2),'r-');
%         plot((-2:11)-0.5,sums{i,3}(:,2,2),'b-');
%         arrayfun(@(x) xline(x,'k:'),[0,1,4,5,7,8]);
%         waitfor(fh)
%     end
% end
% 
% 
% for i=1:size(sums,1)
%     if nnz([sums{i,3}(:,5,1);sums{i,3}(:,5,2)]<0.01)>2
% %         disp(sums{i,3})
% 
%         fh=figure();
%         hold on
% 
%         arrayfun(@(x) xline(x,'k:'),[0,1,4,5,7,8]);
%         waitfor(fh)
%     end
% end




