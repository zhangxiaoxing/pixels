function out=get_stats_by_mem_type(stats,reg,level,opt)
arguments
    stats (1,1) struct {mustBeNonempty}
    reg (1,:) char {mustBeMember(reg,{'between','within','High2Low','Low2High','progres','regres','same_loop','no_loop'})}
    level (1,1) double {mustBePositive,mustBeInteger,mustBeLessThanOrEqual(level,6)}
    opt.congru (1,1) logical = true
    opt.nonmem (1,1) logical = true
    opt.incong (1,1) logical = true
    opt.pvsst (1,1) logical = false
    opt.prog_cross_reg (1,1) logical = false
end
persistent ratiomap OBM1map rings
if isempty(ratiomap) || isempty(OBM1map) || isempty('rings')
    [~,~,ratiomap]=ref.get_pv_sst();
    fstr=load('OBM1map.mat','OBM1map');
    OBM1map=fstr.OBM1map;
    fstr=load(fullfile('bzdata','rings_bz.mat'),'rings');
    rings=fstr.rings;
end
trialtype=char(fieldnames(stats));

% TCOM
if ismember(reg,{'progres','regres'})
    com_str=wave.get_com_map();
    progres=false(size(stats.(trialtype).sess));
    regres=progres;
    for si=reshape(unique(stats.(trialtype).sess),1,[])
        onesess=com_str.(strjoin({'s',num2str(si)},''));
        if opt.prog_cross_reg
            sessidx=find(stats.(trialtype).sess==si & stats.(trialtype).diff_reg(:,level));
        else
            sessidx=find(stats.(trialtype).sess==si);
        end
        for pairidx=reshape(sessidx,1,[])
            %         disp(pairidx)
            if all(stats.(trialtype).mem_type(pairidx,:)==2)
                if all(arrayfun(@(x) onesess.s1.isKey(x), stats.(trialtype).sess_suids(pairidx,:)))
                    if diff(arrayfun(@(x) onesess.s1(x), stats.(trialtype).sess_suids(pairidx,:)),1,2)>0
                        progres(pairidx)=true;
                    else
                        regres(pairidx)=true;
                    end
                end
            elseif all(stats.(trialtype).mem_type(pairidx,:)==4)
                if all(arrayfun(@(x) onesess.s2.isKey(x), stats.(trialtype).sess_suids(pairidx,:)))
                    if diff(arrayfun(@(x) onesess.s2(x), stats.(trialtype).sess_suids(pairidx,:)),1,2)>0
                        progres(pairidx)=true;
                    else
                        regres(pairidx)=true;
                    end
                end
            end
        end
        disp([si,nnz(progres),nnz(regres)]);
    end
end
%loops
if ismember(reg,{'same_loop','no_loop'})
    sameloop=false(size(stats.(trialtype).sess));
    noloop=sameloop;
    for si=reshape(unique(stats.(trialtype).sess),1,[])
        for midx=1:3
            onesess=rings{si,midx};
            if isempty(onesess)
                continue
            end
            sessidx=find(stats.(trialtype).sess==si);
            for pairidx=reshape(sessidx,1,[])
                inloopcnt=sum(onesess==stats.(trialtype).sess_suids(pairidx,1) | onesess==stats.(trialtype).sess_suids(pairidx,2),2);
                if any(inloopcnt==2)
                    sameloop(pairidx)=true;
                else
                    noloop(pairidx)=true;
                end
            end
        end
    end
end



switch reg
    case 'between'
        regsel=false(size(stats.(trialtype).sess));
        for i=1:numel(regsel)
            if opt.pvsst && stats.(trialtype).diff_reg(i,level) && all(cellfun(@(x) ratiomap.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)
                regsel(i)=true;
            elseif ~opt.pvsst && stats.(trialtype).diff_reg(i,level) && all(cellfun(@(x) OBM1map.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)
                regsel(i)=true;
            end
        end
    case 'within'
        regsel=false(size(stats.(trialtype).sess));
        for i=1:numel(regsel)
            if opt.pvsst && stats.(trialtype).same_reg(i,level) && all(cellfun(@(x) ratiomap.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)
                regsel(i)=true;
            elseif ~opt.pvsst && stats.(trialtype).same_reg(i,level) && all(cellfun(@(x) OBM1map.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)
                regsel(i)=true;
            end
        end
    case 'progres'
        regsel=progres;
    case 'regres'
        regsel=regres;
    case 'same_loop'
        regsel=sameloop;
    case 'no_loop'
        regsel=noloop;
        
    case 'High2Low'
        %hierachy
        regsel=false(size(stats.(trialtype).sess));
        for i=1:numel(regsel)
            if opt.pvsst && all(cellfun(@(x) ratiomap.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)...
                    && ratiomap(stats.(trialtype).reg{i,1}{level})<ratiomap(stats.(trialtype).reg{i,2}{level})
                regsel(i)=true;
            elseif ~opt.pvsst && all(cellfun(@(x) OBM1map.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)...
                    && OBM1map(stats.(trialtype).reg{i,1}{level})<OBM1map(stats.(trialtype).reg{i,2}{level})
                regsel(i)=true;
            end
        end
    case 'Low2High'
        %hierachy
        regsel=false(size(stats.(trialtype).sess));
        for i=1:numel(regsel)
            if opt.pvsst && all(cellfun(@(x) ratiomap.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)...
                    && ratiomap(stats.(trialtype).reg{i,1}{level})>ratiomap(stats.(trialtype).reg{i,2}{level})
                regsel(i)=true;
            elseif  ~opt.pvsst && all(cellfun(@(x) OBM1map.isKey(x{level}),stats.(trialtype).reg(i,:)),2)...
                    && all(cellfun(@(x) strcmp(x{2},'CTX'),stats.(trialtype).reg(i,:)),2)...
                    && OBM1map(stats.(trialtype).reg{i,1}{level})>OBM1map(stats.(trialtype).reg{i,2}{level})
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