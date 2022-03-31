function [is_diff,is_same,h2l,l2h]=diff_at_level(reg,opt)
arguments
    reg
    opt.hierarchy (1,1)  logical = false
    opt.hiermap (1,:) char {mustBeMember(opt.hiermap,{'pvsst','OBM1','AON','LAT'})} = 'AON'
    opt.mincount (1,1) double = 100
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
end
assert(opt.hierarchy)
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');

switch opt.hiermap
    case 'pvsst'
        [~,~,hiermap]=ref.get_pv_sst();
    case 'OBM1'
        fstr=load('OBM1map.mat','OBM1map');
        hiermap=fstr.OBM1map;
    case 'AON'
        anov=sink_src_mat(:,src_ccfid==idmap.reg2ccfid('AON'));
        hiermap=containers.Map(cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
            anov);
    case 'LAT'
        anov=sink_src_mat(:,src_ccfid==idmap.reg2ccfid('LAT'));
        hiermap=containers.Map(cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
            anov);
end

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
if opt.hierarchy
    is_diff=[];
    is_same=false(size(reg,1),6);
    h2l=false(size(reg,1),6);
    l2h=false(size(reg,1),6);
    grey_ccfids=int32(cell2mat(idmap.reg2ccfid.values(ephys.getGreyRegs('mincount',opt.mincount,'range',opt.range))));
    greysel=all(ismember(reg(:,5,:),grey_ccfids),3);
    is_same(:,5)=greysel & reg(:,5,1)==reg(:,5,2);
    idmap.ccfid2reg(0)={'MISSING'};
    reg_cell=cellfun(@(x) x,idmap.ccfid2reg.values(squeeze(num2cell(reg(:,5,:)))));
    rowsel=greysel & all(ismember(reg_cell,hiermap.keys()),2);
    mapv=cell2mat(hiermap.values(reg_cell(rowsel,:)));
    l2h(rowsel,5)=mapv(:,1)<mapv(:,2);
    h2l(rowsel,5)=mapv(:,1)>mapv(:,2);
else
    is_diff=false(size(reg,1),6);
    is_same=false(size(reg,1),6);
    h2l=[];
    l2h=[];
    greysel=all(reg(:,1,:)==567 | reg(:,1,:)==343,3);
    for dep=1:6
        is_diff(:,dep)=greysel & reg(:,dep,1)~=reg(:,dep,2);
        is_same(:,dep)=greysel & reg(:,dep,1)==reg(:,dep,2);
    end
end
end