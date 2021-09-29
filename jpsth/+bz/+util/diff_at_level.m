function [is_diff,is_same,h2l,l2h]=diff_at_level(reg,opt)
arguments
    reg
    opt.hierarchy (1,1)  logical = false
    opt.pvsst (1,1) logical = false
end
persistent ratiomap idmap OBM1map

if isempty(ratiomap) || isempty(idmap) || isempty(OBM1map)
    [~,~,ratiomap]=ref.get_pv_sst();
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
    fstr=load('OBM1map.mat','OBM1map');
    OBM1map=fstr.OBM1map;
end
if opt.hierarchy
    is_diff=[];
    is_same=false(size(reg,1),6);
    h2l=false(size(reg,1),6);
    l2h=false(size(reg,1),6);
    graysel=all(reg(:,2,:)==688,3);
    is_same(:,5)=graysel & reg(:,5,1)==reg(:,5,2) & reg(:,5,1)>0;
    for ri=1:size(reg,1)
        % CTX and labeled at ccf depth 7
        if any(reg(ri,5,:)==0,'all') || any(reg(ri,2,:)~=688,'all') || reg(ri,5,1)==reg(ri,5,2)
            continue
        end
        if ~opt.pvsst && ~all(arrayfun(@(x) OBM1map.isKey(char(idmap.ccfid2reg(x))),squeeze(reg(ri,5,:))))
            disp(arrayfun(@(x) char(idmap.ccfid2reg(x)),squeeze(reg(ri,5,:)),'UniformOutput',false));
            continue
        end
        
        if opt.pvsst
            dhier=diff(arrayfun(@(x) ratiomap(char(idmap.ccfid2reg(x))),squeeze(reg(ri,5,:))));
        else
            dhier=diff(arrayfun(@(x) OBM1map(char(idmap.ccfid2reg(x))),squeeze(reg(ri,5,:))));
        end
        if dhier>0
        	l2h(ri,5)=true;
        else
            h2l(ri,5)=true;
        end
    end

else
    is_diff=false(size(reg,1),6);
    is_same=false(size(reg,1),6);
    h2l=[];
    l2h=[];
    graysel=all(reg(:,1,:)==567 | reg(:,1,:)==343,3);
    for dep=1:6
        is_diff(:,dep)=graysel & reg(:,dep,1)~=reg(:,dep,2);
        is_same(:,dep)=graysel & reg(:,dep,1)==reg(:,dep,2);
    end
end
end