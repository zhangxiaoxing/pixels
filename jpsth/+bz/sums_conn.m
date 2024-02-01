function sums_conn_str=sums_conn(su_meta,opt)
arguments
    su_meta
    opt.poolsize (1,1) double {mustBeInteger,mustBePositive}= 2
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','Naive','any'})} = 'WT'
    opt.inhibit (1,1) logical = false
end
pool = gcp('nocreate');
if isunix && isempty(pool)
    pool=parpool(opt.poolsize);
end
if strcmp(opt.type,'neupix')
    upath=unique(su_meta.allpath);
    fl=fullfile(upath,'bz_corr.mat');
    if strcmp(opt.criteria,'Learning')
        sfn=fullfile('binary','sums_conn_learning.mat');
    elseif strcmp(opt.criteria,'Naive')
        sfn=fullfile('binary','sums_conn_naive.mat');
    else
        if opt.inhibit
            sfn='sums_conn_inhibit.mat';
        else
            sfn=fullfile('binary','sums_conn.mat');
        end
    end
else
    % fl=dir(fullfile('K:','neupix','AIOPTO','BZPART','BZ_XCORR_duo_f*.mat'));
    % sfn='aiopto_sums_conn.mat';
    keyboard();
end
tic
if isunix
    futures=parallel.FevalFuture.empty(numel(fl),0);
    for task_idx = 1:numel(fl)
        futures(task_idx) = parfeval(pool,@sum_one,1,fl(task_idx)); % async significant functional coupling map->reduce
    end
    sums_conn_str=fetchOutputs(futures);
elseif ispc
    out_idx=1;
    for task_idx = 1:numel(fl)
         tmp = sum_one(fl{task_idx}); % async significant functional coupling map->reduce
         if rem(task_idx,10)==0
             disp(task_idx)
         end
         if isfield(tmp,'sig_con') && ~isempty(tmp.sig_con)
             sums_conn_str(out_idx)=tmp;
             out_idx=out_idx+1;
         end
    end
end
toc
blame=vcs.blame();
save(sfn,'sums_conn_str','blame')
end

function out=sum_one(f)
arguments
    f (1,:) char
end
out=struct();
fstr=load(f);
if ~isfield(fstr,'all_probe_ccg')
    return
end
suid=fstr.all_probe_ccg.completeIndex(:,2);
out.sig_con=suid(fstr.all_probe_ccg.sig_con);
out.folder=regexp(f,'.*(?=[\\/]bz_corr.mat)','match','once');
out.ccg_sc=[];
sigccg=cell2mat(arrayfun(@(x) fstr.all_probe_ccg.ccgR(:,fstr.all_probe_ccg.sig_con(x,1),fstr.all_probe_ccg.sig_con(x,2)),...
    1:size(fstr.all_probe_ccg.sig_con,1),'UniformOutput',false));
for i=1:size(fstr.all_probe_ccg.sig_con,1)
    out.qc(i,:)=bz.good_ccg(sigccg(:,i));
    out.ccg_sc=[out.ccg_sc;sigccg(:,i).'];
end
end



