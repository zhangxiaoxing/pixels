load(fullfile('bzdata','rings_bz.mat'),'rings');
cnt=0;
for sid=1:size(rings,1)
    for rsize=3:5
        if isempty(rings{sid,rsize-2}), continue;end
        rmax=size(rings{sid,rsize-2},1);
        n100=ceil(rmax./100);
        for hh=0:(n100-1)
            system (sprintf('bash ring_stats.sh %d %d %d %d',sid,rsize,hh*100+1,min((hh+1)*100,rmax)))
            cnt=cnt+1;
        end
    end
end

function sumsup
    sums_all=cell(1,3);
    for rsize=3:5
       fs=dir(fullfile('bzdata',sprintf('ring_stats_%d*.mat',rsize)));
       sums=cell(0);
       for fi=1:numel(fs)
           fstr=load(fullfile(fs(fi).folder,fs(fi).name),'sums');
           sums=[sums;fstr.sums];
       end
       sums_all{rsize-2}=sums;
    end
    save(fullfile('bzdata','sums_ring_stats_all.mat'),'sums_all');
end