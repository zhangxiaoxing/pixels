pt=load('ring_partial','out');
load(fullfile('bzdata','rings_bz.mat'),'rings');
cnt=0;
for fi=1:size(pt.out,1)
    sid=pt.out{fi,1}(1);
    rsize=pt.out{fi,1}(2);
    sessring=1:size(rings{sid,rsize-2},1);
    if ~isempty(pt.out{fi,2})
        partring=pt.out{fi,2}(1:end-1);
        sessring=sessring(~ismember(sessring,partring));
    end
    for rid=reshape(sessring,1,[])
%         system(sprintf('bash ring_stats.sh %d %d %d',sid,rsize,rid))
        cnt=cnt+1;
    end
end