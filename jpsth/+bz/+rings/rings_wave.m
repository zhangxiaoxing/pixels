% deprecated wave definition
% only evaluate SU-waveid, spike-time independent

% nonmem,overlap,independent,within,cross,ringsize,session,cids,[nan]
function [out,supp_stats]=rings_wave(opt)
arguments
    opt.shufid double {mustBeScalarOrEmpty} = []
    opt.early (1,1) logical = false
    opt.ctx_sel (1,1) logical =true
end
persistent meta rings_shuf
if isempty(meta)
    meta=ephys.util.load_meta();
end
if ~isempty(opt.shufid) && isempty(rings_shuf)
    load(fullfile('bzdata','rings_bz_shuf.mat'),'rings_shuf');
end

if isempty(opt.shufid)
    load(fullfile('bzdata','rings_bz.mat'),'rings'); % data generated with ring_list_bz.m
else
    rings=rings_shuf{opt.shufid};
end

out=[];
supp_stats=[];

for fi=1:size(rings,1)
    for rsidx=1:3
        if isempty(rings{fi,rsidx})
            continue
        end
        sess=fi;
        sesssel=meta.sess==sess;
        cids=meta.allcid(sesssel);
        regs=containers.Map(cids,meta.reg_tree(5,sesssel));
        reg_class=containers.Map(cids,meta.reg_tree(2,sesssel));
        %     memtypes=meta.mem_type(sesssel);
        wave_ids=ephys.get_wave_id(fi,rings{fi,rsidx},'early',opt.early);
        ctx_sel=all(strcmp(reg_class.values(num2cell(rings{fi,rsidx})),'CTX'),2);
        if opt.ctx_sel
            reg_sel=ctx_sel;
        else
            reg_sel=true(size(rings{fi,rsidx},1));
        end

        reg_tagged=~strcmp(regs.values(num2cell(rings{fi,rsidx})),'');
        within_sel=reg_sel & all(reg_tagged,2) & all(reg_tagged(:,2:end)==reg_tagged(:,1),2);
        cross_sel=reg_sel & all(reg_tagged,2) & any(reg_tagged(:,2:end)~=reg_tagged(:,1),2);
        supp_stats=[supp_stats;ctx_sel,...%ctx
            reg_sel & any(wave_ids>0,2) & any(wave_ids==0,2),...%mem-nonmem
            reg_sel & any(ismember(wave_ids,[1 3 5]),2) & any(ismember(wave_ids,[2 4 6]),2) & all(wave_ids~=0,2),... % different sample
            reg_sel & (all(ismember(wave_ids,[1 3 5]),2) | all(ismember(wave_ids,[2 4 6]),2)),...% same sample, if early = false;
            reg_sel & (all(wave_ids==5,2) | all(wave_ids==6,2))];% both
        % nonmem,overlap,independent,within,cross,ringsize,session,cids,[nan]
        out=[out;...
            all(wave_ids==0,2) & reg_sel,...
            reg_sel & all(ismember(wave_ids,[1,5]),2) | all(ismember(wave_ids,[2,6]),2) | all(ismember(wave_ids,[3,5]),2) | all(ismember(wave_ids,[4,6]),2),...
            reg_sel & (any(wave_ids==1,2) & any(wave_ids==3,2) & all(ismember(wave_ids,[1 3]),2)) | (any(wave_ids==2,2) & any(wave_ids==4,2) & all(ismember(wave_ids,[2 4]),2)),...
            within_sel,...
            cross_sel,...
            (rsidx+2)*ones(size(within_sel)),...
            fi*ones(size(within_sel)),...
            rings{fi,rsidx},...
            nan(numel(within_sel),3-rsidx)];
    end
end

