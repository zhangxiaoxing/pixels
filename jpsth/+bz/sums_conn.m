function sums_conn_str=sums_conn(opt)
arguments
    opt.poolsize (1,1) double {mustBeInteger,mustBePositive}= 2
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.prefix (1,:) char = 'BZWT';
    opt.inhibit (1,1) logical = false;
    opt.ccg_strict_qc (1,1) logical = false;
end
pool = gcp('nocreate');
if isempty(pool)
    pool=parpool(opt.poolsize);
end
if strcmp(opt.type,'neupix')
    if strcmp(opt.criteria,'Learning')
        fl=dir(fullfile('bzdata',sprintf('%s_BZ_XCORR_duo_f*_Learning.mat',opt.prefix)));
        sfn='sums_conn_learning.mat';
    else
        if ispc
            fl=dir(fullfile('bzdata',sprintf('%s_BZ_XCORR_duo_f*.mat',opt.prefix)));
        elseif isunix
            fl=dir(fullfile('/media/SSD2','bzdata',sprintf('%s_BZ_XCORR_duo_f*.mat',opt.prefix)));
        end
        if opt.inhibit
            sfn='sums_conn_inhibit.mat';
        else
            sfn='sums_conn.mat';
        end
    end
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
out.ccg_sc=[];
sigccg=cell2mat(arrayfun(@(x) fstr.mono.ccgR(:,fstr.mono.sig_con(x,1),fstr.mono.sig_con(x,2)),...
    1:size(fstr.mono.sig_con,1),'UniformOutput',false));
for i=1:size(fstr.mono.sig_con,1)
    out.qc(i,:)=bz.good_ccg(sigccg(:,i));
    if out.qc(i,1)>0 ... % peak dir
            && out.qc(i,3)<10 ... % noise
            && out.qc(i,4)<5 ... % fwhm
            && out.qc(i,2)>253 ... % peak time
            && out.qc(i,2)<257 % peak time
        out.ccg_sc=[out.ccg_sc;...
            fstr.mono.sig_con(i,:),out.qc(i,:),sigccg(:,i)];
    end
end
end



