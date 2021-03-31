% pool=parpool(2);
pool=gcp(); % assuming available pool
prefix='0315';
fl=dir(sprintf('%s_BZ_XCORR_duo_f*.mat',prefix));
futures=parallel.FevalFuture.empty(numel(fl),0);

for task_idx = 1:numel(fl)
    futures(task_idx) = parfeval(pool,@sum_one,1,fl(task_idx)); % async significant functional coupling map->reduce
end
sums_conn_str=fetchOutputs(futures);
save('sums_conn.mat','sums_conn_str')


function out=sum_one(f)
fstr=load(fullfile(f.folder,f.name));
suid=fstr.mono.completeIndex(:,2);
out.sig_con=suid(fstr.mono.sig_con);
out.folder=fstr.folder;
% out.ccg=cell2mat(arrayfun(@(x) fstr.mono.ccgR(:,fstr.mono.sig_con(x,1),fstr.mono.sig_con(x,2)),1:size(fstr.mono.sig_con,1),'UniformOutput',false));
end