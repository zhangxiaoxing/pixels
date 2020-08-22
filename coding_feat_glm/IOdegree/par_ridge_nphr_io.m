function results=par_ridge_nphr_selectivity(effs, tasks,lambda,use_mean,early_late)
addpath('FastRidge');

if ~exist('lambda','var')
   lambda=0.1; 
end
if ~exist('use_mean','var')
    use_mean=true;
end

transFrac=transientFraction();%early col3, late col5
iostats=ioDegree(); %bin 1-6
% iodiff=ioDiff();
dec_stats=dec_accuracy();%entire col 1, early col 2, late col3
opgenStats=hemOpgen(use_mean,'nphr');
results=cell(size(tasks));
for t=tasks
    results{t}=one_corr(effs,t, lambda, use_mean, transFrac,iostats,opgenStats,dec_stats,early_late);
end
return
% p = gcp();
% futures=parallel.FevalFuture.empty(0);
% for use_eff=effs
%     for task_idx=tasks
%         futures(end+1) = parfeval(p,@one_corr,0,use_eff,task_idx, lambda, use_mean, subset,transFrac,iostats); % Square size determined by idx
%     end
% end
% % Collect the results as they become available.
% for idx = 1:numel(futures)
%     % fetchNext blocks until next results are available.
%     completedIdx = fetchNext(futures);
%     fprintf('Got result with index: %d.\n', completedIdx);
% end

end


function int_result=one_corr(use_eff,task_idx,lambda,use_mean,...
    transFrac,iostats,opgenStats,dec_stats,early_late)

glm_mat=[];
regions=cell(0,2);

for i=1:size(opgenStats.Mapping,1)
    opgen_reg=opgenStats.Mapping{i,1};
    ephys_reg=opgen_reg;
    opgen_idx=i;
    frac_idx=strcmp(ephys_reg,transFrac.reg);
    io_idx=strcmp(ephys_reg,iostats.reg);
    dec_idx=strcmp(ephys_reg,dec_stats.reg);
    if nnz(opgen_idx) && nnz(frac_idx) && nnz(io_idx) && nnz(dec_idx)...
            && (transFrac.stats(frac_idx,1)>=50)
%         if ~nnz(iodiff_idx)
%             disp(ephys_reg);
%             keyboard
%             
%         end
        
        opgen=opgenStats.value_arr(opgen_idx,task_idx);
        if ~exist('early_late','var') || strcmp(early_late,'entire')
            disp("not ready")
            keyboard
            frac=transFrac.stats(frac_idx,[3 5]);% col 3 for early col 5 for late
            iodeg=[mean(iostats.dens(io_idx,1:3)),mean(iostats.dens(io_idx,4:6))]; % % 6 bins 
            dec_diff=dec_stats.stats(dec_idx,[5 6]); % 5 for early decoding 6 for late
%             iodiff_val=iodiff.stats(iodiff_idx,1);
        elseif strcmp(early_late,'early')
            frac=transFrac.stats(frac_idx,3);% col 3 for early col 5 for late
            iodeg=iostats.dens(io_idx,[4:6,11,14]); % % 1-9 (in,out,local) in (entire early late),10-12, io diff, 13-15, io-diff-idx
            dec_diff=dec_stats.stats(dec_idx,[2 5]);
%             iodiff_val=iodiff.stats(iodiff_idx,2);
        elseif strcmp(early_late,'late')
            frac=transFrac.stats(frac_idx,5);% col 3 for early col 5 for late
            iodeg=iostats.dens(io_idx,[7:9,12,15]); % % 6 bins 
            dec_diff=dec_stats.stats(dec_idx,[3 6]); % 5 for early decoding 6 for late
%             iodiff_val=iodiff.stats(iodiff_idx,3);
        else
            disp('error early late delay parameter');
            keyboard
        end
        if any(isnan([opgen,frac,iodeg,dec_diff]))
            continue
        end

        glm_mat(end+1,:)=[opgen,frac,dec_diff,iodeg];
        regions(end+1,:)={opgen_reg,ephys_reg};
    elseif transFrac.stats(frac_idx,1)>=50
%         disp(ephys_reg)
%         disp([nnz(opgen_idx),nnz(frac_idx), nnz(io_idx),nnz(dec_idx),nnz(iodiff_idx)]);
%         keyboard;
    end

end

disp(size(glm_mat))

int_result=cell(0,7);

for n=1:(size(glm_mat,2)-1)
    C=nchoosek(2:size(glm_mat,2),n);
    for i=1:size(C,1)
        [beta, b0, tau2, DOF, lambda, score] = fastridge(glm_mat(:,C(i,:)), glm_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
%         mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'linear');
        [rr,pp]=corr(glm_mat(:,1),glm_mat(:,C(i,:))*beta,'type','Spearman');
        new_idx=size(int_result,1)+1;
        int_result{new_idx,1}=beta;
        int_result{new_idx,2}=score;
%         int_result{new_idx,3}=rr(1,2)^2;
        int_result{new_idx,3}=rr;
%         int_result{new_idx,4}=pp(1,2);
        int_result{new_idx,4}=pp;
        int_result{new_idx,5}='linear';
        int_result{new_idx,6}=C(i,:);
        int_result{new_idx,7}=[glm_mat(:,1),glm_mat(:,C(i,:))*beta];
    end
end
if false
    for n=2:4
        C=nchoosek(2:size(glm_mat,2),n);
        for i=1:size(C,1)
            X_from=glm_mat(:,C(i,:));
            X=x2fx(X_from,'interaction');
            X(:,1)=[];
            [beta, b0, tau2, DOF, lambda, score] = fastridge(X, glm_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
    %         mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'linear');
            [rr,pp]=corrcoef(glm_mat(:,1),X*beta);
            new_idx=size(int_result,1)+1;
            int_result{new_idx,1}=beta;
            int_result{new_idx,2}=score;
            int_result{new_idx,3}=rr(1,2)^2;
            int_result{new_idx,4}=pp(1,2);
            int_result{new_idx,5}='interact';
            int_result{new_idx,6}=C(i,:);
            int_result{new_idx,7}=[glm_mat(:,1),X*beta];
        end
    end
end
%% cross validation
% [~,low_aic_idx]=min(cell2mat(int_result(:,2)));
% cv_results=nan(size(glm_mat,1),2);
% for i=1:size(glm_mat,1)
%     cv_mat=glm_mat;
%     cv_mat(i,:)=[];
%     if strcmp(int_result{low_aic_idx,5},'linear')
%         X=cv_mat(:,int_result{low_aic_idx,6});
%         [beta, b0, tau2, DOF, lambda, score] = fastridge(X, cv_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
%         X=glm_mat(i,int_result{low_aic_idx,6});
%         pred=X*beta+b0;
%         cv_results(i,:)=[glm_mat(i,1),pred]; % target prediction
%     else
%         X_from=cv_mat(:,int_result{low_aic_idx,6});
%         X=x2fx(X_from,'interaction');
%         X(:,1)=[];
%         [beta, b0, tau2, DOF, lambda, score] = fastridge(X, cv_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
%         X_from=glm_mat(i,int_result{low_aic_idx,6});
%         X=x2fx(X_from,'interaction');
%         pred=X(2:end)*beta+b0;
%         cv_results(i,:)=[glm_mat(i,1),pred]; % target prediction
%     end
% end
% figure()
% plot(cv_results(:,1),cv_results(:,2),'k.')
% ylim([-1,1])
% xlim([-1,1])
% [r,p]=corrcoef(cv_results(:,1),cv_results(:,2));
% r=r(1,2);
% p=p(1,2);
if true
    value_labels={'early correct','early hit','early CR','late correct','late hit','late CR',};

    suffix=value_labels{task_idx};
    if use_eff
        suffix=['EFF_',suffix];
    end
    if use_mean
        suffix=['Mean_',suffix];
    else
        suffix=['Median_',suffix];
    end
    save(sprintf('Ridge_frac_iodeg_%s_lambda_%0.2f_%s_%s.mat',opgenStats.opsin,lambda,suffix,early_late),'int_result','lambda','regions')
end
end

%%
function transFrac=transientFraction()
datapath='k:\code\transient_6.hdf5';
%sust, transient, switched, unclassified, early_in_6s, late_in_6s, prefer_s
sus_trans = h5read(datapath,'/sus_trans');
reg = h5read(datapath,'/reg');
reg=cellfun(@(x) deblank(x),reg,'UniformOutput',false);
reg_set=unique(reg);
reg_set=reg_set(~strcmp(reg_set,'Unlabeled') & ~strcmp(reg_set,'root'));


stats=[];
for one_reg=reg_set'
    reg_sel=strcmp(reg,strtrim(one_reg));
%     keyboard
    count=nnz(reg_sel);
%     trans=nnz(reg_sel & sus_trans(:,2));
%     sust=nnz(reg_sel & sus_trans(:,1));
    early=nnz(reg_sel & sus_trans(:,5));
    late=nnz(reg_sel & sus_trans(:,6));
    stats(end+1,:)=([count,early,early/count,late,late/count]);
end
transFrac=struct();
transFrac.reg=reg_set;
transFrac.stats=stats;
% dist=[x[3] for x in sums if x[1]>=100]
% dmin=np.floor(np.min(dist)*100).astype(np.int)
% dmax=np.floor(np.max(dist)*100).astype(np.int)
end

%%
function iostats=ioDegree()
iostats=struct();
% load('k:\code\jpsth\reg_keep.mat')
% dens=nan(140,6);
% for bin=1:6
% load(sprintf('k:\\code\\jpsth\\0626_selec_conn_mat_duo_6s_%d_%d.mat',bin,bin+1))
% load(sprintf('k:\\code\\jpsth\\0810_selec_pair_mat_duo_6s_%d_%d.mat',bin,bin+1))
% for i=1:length(reg_set)
%     outPair=sum(pair_mat(:,i));
%     outConn=sum(conn_mat_S1(:,i));
%     if outPair>0
%         dens(i,bin)=outConn/outPair;
%     end
% end
% end
load('k:\code\jpsth\io_sel.mat');
reg_sel=io_early_delay(:,1)>=100 & io_late_delay(:,1)>=100 & io_entire_delay(:,1)>=100 & ...
    io_early_delay(:,10)>=100 & io_late_delay(:,10)>=100 & io_entire_delay(:,10)>=100;

iostats.reg=reg_set(reg_sel);
stats=[io_entire_delay(reg_sel,[3 7 12]),io_early_delay(reg_sel,[3 7 12]),io_late_delay(reg_sel,[3 7 12])];
iodiff=[diff(io_entire_delay(reg_sel,[3,7]),1,2),diff(io_early_delay(reg_sel,[3,7]),1,2),diff(io_late_delay(reg_sel,[3,7]),1,2)];
iodfIdx=[iodiff(:,1)./sum(io_entire_delay(reg_sel,[3,7]),2),...
    iodiff(:,1)./sum(io_early_delay(reg_sel,[3,7]),2),...
    iodiff(:,1)./sum(io_late_delay(reg_sel,[3,7]),2)];
iostats.dens=[stats,iodiff,iodfIdx];
% keyboard
iostats.reg{strcmp(iostats.reg,'SSp')}='SSp-bfd';
iostats.reg{strcmp(iostats.reg,'DG')}='DG-mo';
end

% function iodiff=ioDiff()
% load('k:\code\jpsth\io_sel.mat'); %save('io_sel.mat','ioselstats','io_entire_delay','io_early_delay','io_late_delay','reg_set');
% %     in_out_sel(reg_idx,:)=[pair_count,in_conn_S1,in_conn_S1/pair_count, ...%1 2 3
% %         in_sel_S1,in_sel_S1/pair_count,...% 4 5
% %         out_conn_S1,out_conn_S1/pair_count,...% 6 7
% %         out_sel_S1,out_sel_S1/pair_count,...% 8 9
% %         auto_pair,auto_conn_S1,auto_conn_S1/auto_pair]; % 10 11 12
% iodiff=struct();
% iodiff.reg=reg_set;
% iodiff.stats=[diff(io_entire_delay(:,[3,7]),1,2),diff(io_early_delay(:,[3,7]),1,2),diff(io_late_delay(:,[3,7]),1,2)];
% iodiff.reg{strcmp(iodiff.reg,'SSp')}='SSp-bfd';
% iodiff.reg{strcmp(iodiff.reg,'DG')}='DG-mo';
% end


%%
function dec_stats=dec_accuracy()
dec_stats=struct();
fpath='K:\code\jpsth';
%ref:K:\code\jpsth\io_sel_ctd_accuracy.m
% delay corr, early corr, late corr, delay err, early err, late err
load(fullfile(fpath,'ctd_accuracy.mat'));
dec_stats.reg=c_reg;
%entire delay, early, late
dec_stats.stats=[dec_accu(:,1)-dec_accu(:,4),dec_accu(:,2)-dec_accu(:,5),dec_accu(:,3)-dec_accu(:,6),...
    dec_accu(:,1),dec_accu(:,2),dec_accu(:,3)];
end

%%
function opgenStats=hemOpgen(use_mean,use_vgat)
if exist('use_vgat','var') && strcmp(use_vgat,'vgat')
     perfpath='K:\code\coding_feat_glm\IOdegree\vgat_new\optogenetic';
     load(fullfile(perfpath,'anaData_ED.mat'));
     value_arrE=cell2mat(cellfun(@(x) mean(x(:,1)),effectSize,'UniformOutput',false));
     load(fullfile(perfpath,'anaData_LD.mat'));
     value_arrL=cell2mat(cellfun(@(x) mean(x(:,1)),effectSize,'UniformOutput',false));
     opgenStats.value_arr=[value_arrE,value_arrL];
     opgenStats=struct();
     opgenStats.value_arr=[value_arrE,value_arrL];
     opgenStats.value_arr=opgenStats.value_arr(3:end,:);
     opgenStats.Mapping=testID;
     opgenStats.opsin='vgat';
else
    %% 0711 update import perf from HEM result
    perfpath='K:\code\coding_feat_glm\IOdegree\ZX\mapping NpHR\';
    load(fullfile(perfpath,'effectsize_ED.mat'));
    if use_mean
        value_arrE=cell2mat(cellfun(@(x) mean(x), effectsize_normalized,'UniformOutput',false));
    else
        value_arrE=cell2mat(cellfun(@(x) median(x), effectsize_normalized,'UniformOutput',false));
    end
    load(fullfile(perfpath,'effectsize_LD.mat'));
    if use_mean
        value_arrL=cell2mat(cellfun(@(x) mean(x), effectsize_normalized,'UniformOutput',false));
    else
        value_arrL=cell2mat(cellfun(@(x) median(x), effectsize_normalized,'UniformOutput',false));
    end
    %%

    opgenStats=struct();
    opgenStats.value_arr=[value_arrE,value_arrL];

%     if false
%         scatter(value_arrE(:,1),value_arrL(:,1))
%         [r,p]=corr(value_arrE(:,1),value_arrL(:,1),'type','Spearman')
%     end
    opgenStats.Mapping=Mapping;
    opgenStats.opsin='nphr';
end
end
