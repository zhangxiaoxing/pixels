function out=get_stats_by_mem_type(stats,reg,level,opt)
arguments
    stats (1,1) struct {mustBeNonempty}
    reg (1,:) char {mustBeMember(reg,{'between','within'})}
    level (1,1) double {mustBePositive,mustBeInteger,mustBeLessThanOrEqual(level,6)}
    opt.congru (1,1) logical = true
    opt.nonmem (1,1) logical = true
    opt.incong (1,1) logical = true
end
trialtype=char(fieldnames(stats));
if strcmp(reg,'between')
    regsel=stats.(trialtype).diff_reg(:,level);
else
    regsel=stats.(trialtype).same_reg(:,level);
end

if opt.congru
    out.congru=stats.(trialtype).coeff(...
        (all(ismember(stats.(trialtype).mem_type,1:2),2) | all(ismember(stats.(trialtype).mem_type,3:4),2))...
        & regsel & ~stats.(trialtype).skip & stats.(trialtype).rsq>0,:);
end
if opt.nonmem
    out.nonmem=stats.(trialtype).coeff(...
        all(stats.(trialtype).mem_type==0,2)...
        & regsel & ~stats.(trialtype).skip & stats.(trialtype).rsq>0,:);
end

if opt.incong
    out.incong=stats.(trialtype).coeff(...
        any(ismember(stats.(trialtype).mem_type,1:2),2)...
        & any(ismember(stats.(trialtype).mem_type,3:4),2)...
        & regsel & ~stats.(trialtype).skip & stats.(trialtype).rsq>0,:);
end


end