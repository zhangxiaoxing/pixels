function sums_conn_str=sums_conn(opt)
arguments
    opt.poolsize (1,1) double {mustBeInteger,mustBePositive}= 2
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.prefix (1,:) char = '0315';
end
pool = gcp('nocreate');
if isempty(pool)
    pool=parpool(opt.poolsize);
end
if strcmp(opt.type,'neupix')
    fl=dir(sprintf('%s_BZ_XCORR_duo_f*.mat',opt.prefix));
    sfn='sums_conn.mat';
else
    fl=dir(fullfile('K:','neupix','AIOPTO','BZPART','BZ_XCORR_duo_f*.mat'));
    sfn='aiopto_sums_conn.mat';
end
futures=parallel.FevalFuture.empty(numel(fl),0);

for task_idx = 1:numel(fl)
    futures(task_idx) = parfeval(pool,@sum_one,1,fl(task_idx)); % async significant functional coupling map->reduce
end
sums_conn_str=fetchOutputs(futures);

save(sfn,'sums_conn_str')
end

function out=sum_one(f)
fstr=load(fullfile(f.folder,f.name));
suid=fstr.mono.completeIndex(:,2);
out.sig_con=suid(fstr.mono.sig_con);
out.folder=fstr.folder;
% out.ccg=cell2mat(arrayfun(@(x) fstr.mono.ccgR(:,fstr.mono.sig_con(x,1),fstr.mono.sig_con(x,2)),1:size(fstr.mono.sig_con,1),'UniformOutput',false));
end