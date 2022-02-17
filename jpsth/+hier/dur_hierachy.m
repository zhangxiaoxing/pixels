% preferably load disk data
% dur=wave.duration_wave('ctx',false,'extra_ITI',true);

%% temporary work-around
load('duration_wave.mat','out')
dur=out;
%%

meta=ephys.util.load_meta();
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));

if false % one-time consistency check
    assert(isequal(meta.sess,dur.meta(:,1)));
    assert(isequal(meta.allcid,dur.meta(:,2)));
end

ctxsel=strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),'');
cnusel=strcmp(meta.reg_tree(2,:),'CNU');

% ctxregs=unique(meta.reg_tree(5,ctxsel));
chregs=unique(meta.reg_tree(5,ctxsel|cnusel));
dur_reg_map=containers.Map();
for r=chregs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    pos=nnz(strcmp(meta.reg_tree(5,:),r{1}) & any((dur.wrsp(:,4:end)<0.05).'));
    dur_reg_map(r{1})=[pos/cnt,pos,cnt];
end
% load('OBM1Map.mat','OBM1map');

sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/CH_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/CH_srcs');
src_sink_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');

src_idx_map=containers.Map(src_ccfid,1:numel(src_ccfid)); %{ccfid:proj_dense_mat_idx}
sink_idx_map=containers.Map(sink_ccfid,1:numel(sink_ccfid)); %{ccfid:proj_dense_mat_idx}

if false
    figure();
    imagesc(src_sink_mat);
end

%% efferent_proj_dense(src) ~ dur_proportion

allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
intersect_regs=intersect(allen_regs,chregs);
idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
dur_prop_cell=dur_reg_map.values(intersect_regs);
dur_prop=cellfun(@(x) x(1),dur_prop_cell);

for ii=1:numel(sink_ccfid)
    allen_density=log10(src_sink_mat(ii,idx4corr)); % from idx4corr, to one alternating target
    [r,p]=corr(allen_density.',dur_prop);
    if p<0.05
        figure('Color','w')
        scatter(allen_density,dur_prop,4,'red','filled','o')
        text(allen_density,dur_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
        title(['Allen connectivity from each region to', idmap.ccfid2full(sink_ccfid(ii))])
        xlabel('Projection density (log10)')
        ylabel('Proportion of duration neuron')
        text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    end
end

%% afferent_proj_dense(sink) ~ dur_proportion

allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(sink_ccfid)),'UniformOutput',false);
intersect_regs=intersect(allen_regs,chregs);
idx4corr=cell2mat(sink_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
dur_prop_cell=dur_reg_map.values(intersect_regs);
dur_prop=cellfun(@(x) x(1),dur_prop_cell);

for ii=1:numel(src_ccfid)
    allen_density=log10(src_sink_mat(idx4corr,ii)); % from one alternating target, to idx4corr
    [r,p]=corr(allen_density,dur_prop);
    if p<0.05
        figure('Color','w')
        scatter(allen_density,dur_prop,4,'red','filled','o')
        text(allen_density,dur_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
        title(['Allen connectivity to each region from', idmap.ccfid2full(src_ccfid(ii))])
        xlabel('Projection density (log10)')
        ylabel('Proportion of duration neuron')
        text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    end
end