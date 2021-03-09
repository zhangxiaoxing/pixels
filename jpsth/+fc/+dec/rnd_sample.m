function out=rnd_sample(sums,sel_S1,sel_S2,rpts,opt)
arguments
    sums (:,3) cell
    sel_S1 (:,1) logical
    sel_S2 (:,1) logical
    rpts (1,1) double = 2
    opt.featIdx (1,1) double = 3
    opt.trials (1,1) double = 20
end

binn=size(sums{1,3},3);

s1=find(sel_S1);
s2=find(sel_S2);
out=struct();
for rpt=1:rpts
    s61s=randsample(s1, opt.trials);
    s62s=randsample(s2, opt.trials);
    out.s1(:,:,:,rpt)=cell2mat(arrayfun(@(x) sums{x,3}(opt.featIdx,s61s,:),1:size(sums,1),'UniformOutput',false)');
    out.s2(:,:,:,rpt)=cell2mat(arrayfun(@(x) sums{x,3}(opt.featIdx,s62s,:),1:size(sums,1),'UniformOutput',false)');
end
% keyboard();
