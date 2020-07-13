function results=par_ridge_nphr_selectivity(effs, tasks,lambda,use_mean,subset)
addpath('FastRidge');
if ~exist('subset','var')
    subset=false;
end
if ~exist('lambda','var')
   lambda=0.1; 
end
if ~exist('use_mean','var')
    use_mean=true;
end

transFrac=transientFraction();
iostats=ioDegree();
opgenStats=hemOpgen(use_mean);
results=cell(1,6);
for t=tasks
    results{t}=one_corr(effs,t, lambda, use_mean, subset,transFrac,iostats,opgenStats);
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


function int_result=one_corr(use_eff,task_idx,lambda,use_mean,subset,transFrac,iostats,opgenStats)

glm_mat=[];
regions=cell(0,2);

for i=1:size(opgenStats.Mapping,1)
    opgen_reg=opgenStats.Mapping{i,1};
    ephys_reg=opgen_reg;
    opgen_idx=i;
    frac_idx=strcmp(ephys_reg,transFrac.reg);
    io_idx=strcmp(ephys_reg,iostats.reg);
    if nnz(opgen_idx) && nnz(frac_idx) && nnz(io_idx) && (transFrac.stats(frac_idx,1)>=50)

        opgen=opgenStats.value_arr(opgen_idx,task_idx);
        frac=transFrac.stats(frac_idx,[3 5]);% col 3 for transient col 5 for sust
        iodeg=iostats.dens(io_idx,[2 4]); % % 6 bins 
        if any(isnan([opgen,frac,iodeg]))
            continue
        end
        glm_mat(end+1,:)=[opgen,frac,iodeg];
        regions(end+1,:)={opgen_reg,ephys_reg};
    end

end

disp(size(glm_mat))

int_result=cell(0,7);

for n=1:(size(glm_mat,2)-1)
    C=nchoosek(2:size(glm_mat,2),n);
    for i=1:size(C,1)
        [beta, b0, tau2, DOF, lambda, score] = fastridge(glm_mat(:,C(i,:)), glm_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
%         mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'linear');
        [rr,pp]=corrcoef(glm_mat(:,1),glm_mat(:,C(i,:))*beta);
        new_idx=size(int_result,1)+1;
        int_result{new_idx,1}=beta;
        int_result{new_idx,2}=score;
        int_result{new_idx,3}=rr(1,2)^2;
        int_result{new_idx,4}=pp(1,2);
        int_result{new_idx,5}='linear';
        int_result{new_idx,6}=C(i,:);
        int_result{new_idx,7}=[glm_mat(:,1),glm_mat(:,C(i,:))*beta];
    end
end
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
    save(sprintf('Ridge_frac_iodeg_nphr_lambda_%0.2f_%s.mat',lambda,suffix),'int_result','lambda','regions')
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
load('k:\code\jpsth\reg_keep.mat')
dens=nan(140,6);
for bin=1:6
load(sprintf('k:\\code\\jpsth\\0626_selec_conn_mat_duo_6s_%d_%d.mat',bin,bin+1))
load(sprintf('k:\\code\\jpsth\\0626_selec_pair_mat_duo_6s_%d_%d.mat',bin,bin+1))
for i=1:length(reg_set)
    outPair=sum(pair_mat(:,i));
    outConn=sum(conn_mat_S1(:,i));
    if outPair>0
        dens(i,bin)=outConn/outPair;
    end
end
end
iostats.reg=reg_set;
iostats.dens=dens;
% keyboard
end
%%
function opgenStats=hemOpgen(use_mean)

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
opgenStats.Mapping=Mapping;
end
