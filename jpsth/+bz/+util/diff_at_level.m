function [is_diff,is_same,h2l,l2h,hierv]=diff_at_level(reg,opt)
arguments
    reg
    opt.hierarchy (1,1)  logical = false
    opt.hiermap (1,:) char
    opt.mincount (1,1) double = 100
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    opt.descend (1,1) logical = false
end



if opt.hierarchy
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
    src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
    sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');
    src_reg=cellfun(@(x) x, idmap.ccfid2reg.values(num2cell(src_ccfid)));
    assert(ismember(opt.hiermap,[src_reg;{'pvsst'};{'OBM1'}]),'Unknown hierachical mapping');


    switch opt.hiermap
        case 'pvsst'
            [~,~,hiermap]=ref.get_pv_sst();
            assert(~opt.descend,'desending order not supported');
        case 'OBM1'
            fstr=load('OBM1map.mat','OBM1map');
            hiermap=fstr.OBM1map;
            assert(~opt.descend,'desending order not supported');
        case 'pct'
            hiermap=opt.hiermap_map;
        otherwise
            anov=sink_src_mat(:,src_ccfid==idmap.reg2ccfid(opt.hiermap));
            if opt.descend
                anov=-anov;
            end
            hiermap=containers.Map(cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
                anov);

    end

    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));

    is_diff=[];
    [is_same,h2l,l2h]=deal(false(size(reg,1),6));
%     h2l=false(size(reg,1),6);
%     l2h=false(size(reg,1),6);
    hierv=zeros(size(reg,1),2);
    grey_ccfids=int32(cell2mat(idmap.reg2ccfid.values(ephys.getGreyRegs('mincount',opt.mincount,'range',opt.range))));
    greysel=all(ismember(reg(:,5,:),grey_ccfids),3);
    is_same(:,5)=greysel & reg(:,5,1)==reg(:,5,2);
    idmap.ccfid2reg(0)={'MISSING'};
    reg_cell=cellfun(@(x) x,idmap.ccfid2reg.values(squeeze(num2cell(reg(:,5,:)))));
    rowsel=greysel & all(ismember(reg_cell,hiermap.keys()),2);
    hierv(rowsel,:)=cell2mat(hiermap.values(reg_cell(rowsel,:)));
    l2h(rowsel,5)=hierv(rowsel,1)<hierv(rowsel,2);
    h2l(rowsel,5)=hierv(rowsel,1)>hierv(rowsel,2);

else
    is_diff=false(size(reg,1),6);
    is_same=false(size(reg,1),6);
    [h2l,l2h,hierv]=deal([]);
    switch opt.range
        case 'grey'
            greysel=all(reg(:,1,:)==567 | reg(:,1,:)==343,3);
        case 'CTX'
            greysel=all(reg(:,2,:)==688,3);
        case 'CH'
            greysel=all(reg(:,1,:)==567,3);
    end
    for dep=1:6
        is_diff(:,dep)=greysel & reg(:,dep,1)~=reg(:,dep,2) & all(reg(:,dep,:)>0,3);
        is_same(:,dep)=greysel & reg(:,dep,1)==reg(:,dep,2) & all(reg(:,dep,:)>0,3);
    end
end
end