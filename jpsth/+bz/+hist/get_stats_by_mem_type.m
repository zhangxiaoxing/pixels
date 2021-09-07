function out=get_stats_by_mem_type(stats,reg,level,opt)
arguments
    stats (1,1) struct {mustBeNonempty}
    reg (1,:) char {mustBeMember(reg,{'between','within','High2Low','Low2High'})}
    level (1,1) double {mustBePositive,mustBeInteger,mustBeLessThanOrEqual(level,6)}
    opt.congru (1,1) logical = true
    opt.nonmem (1,1) logical = true
    opt.incong (1,1) logical = true
end
persistent ratiomap
if isempty(ratiomap)
    [~,~,ratiomap]=ref.get_pv_sst();
end
trialtype=char(fieldnames(stats));
switch reg
    case 'between' 
        regsel=stats.(trialtype).diff_reg(:,level);
    case 'within'
        regsel=stats.(trialtype).same_reg(:,level);
    case 'High2Low'
        %hierachy
        regsel=false(size(stats.(trialtype).sess));
        for i=1:numel(regsel)
            if all(cellfun(@(x) ratiomap.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)...
                    && ratiomap(stats.(trialtype).reg{i,1}{level})<ratiomap(stats.(trialtype).reg{i,2}{level})
                regsel(i)=true;
            end
        end
    case 'Low2High'
        %hierachy
        regsel=false(size(stats.(trialtype).sess));
        for i=1:numel(regsel)
            if all(cellfun(@(x) ratiomap.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)...
                    && ratiomap(stats.(trialtype).reg{i,1}{level})>ratiomap(stats.(trialtype).reg{i,2}{level})
                regsel(i)=true;
            end
        end
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