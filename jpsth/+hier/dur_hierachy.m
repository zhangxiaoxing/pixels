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

%% duration
dur_reg_map=containers.Map();
for r=chregs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    pos=nnz(strcmp(meta.reg_tree(5,:),r{1}) & any((dur.wrsp(:,4:end)<0.05).'));
    dur_reg_map(r{1})=[pos/cnt,pos,cnt];
end

%% fraction
waveid=ephys.get_wave_id(meta.sess,meta.allcid);
sens_reg_map=containers.Map();
for r=chregs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    pos=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & waveid>0);
    sens_reg_map(r{1})=[pos/cnt,pos,cnt];
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

feat_reg_map=sens_reg_map;
% feat_reg_map=dur_reg_map;

min_max_list=[];

%% efferent_proj_dense(src) ~ feature_proportion

allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
intersect_regs=intersect(allen_regs,chregs);
idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
feat_prop_cell=feat_reg_map.values(intersect_regs);
feat_prop=cellfun(@(x) x(1),feat_prop_cell);

for ii=1:numel(sink_ccfid)
    allen_density=log10(src_sink_mat(ii,idx4corr)); % from idx4corr, to one alternating target
    [r,p]=corr(allen_density.',feat_prop);
    if p<0.05 && false
        figure('Color','w')
        scatter(allen_density,feat_prop,4,'red','filled','o')
        text(allen_density,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
        title(['Allen connectivity from each region to', idmap.ccfid2full(sink_ccfid(ii))])
        xlabel('Projection density (log10)')
        ylabel('Proportion of coding neuron')
        text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    end
    min_max_list=[min_max_list;1,ii,r,p];
end

%% afferent_proj_dense(sink) ~ feature_proportion

allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(sink_ccfid)),'UniformOutput',false);
intersect_regs=intersect(allen_regs,chregs);
idx4corr=cell2mat(sink_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
feat_prop_cell=feat_reg_map.values(intersect_regs);
feat_prop=cellfun(@(x) x(1),feat_prop_cell);

for ii=1:numel(src_ccfid)
    allen_density=log10(src_sink_mat(idx4corr,ii)); % from one alternating target, to idx4corr
    [r,p]=corr(allen_density,feat_prop);
    if p<0.05 && false
        figure('Color','w')
        scatter(allen_density,feat_prop,4,'red','filled','o')
        text(allen_density,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
        title(['Allen connectivity to each region from', idmap.ccfid2full(src_ccfid(ii))])
        xlabel('Projection density (log10)')
        ylabel('Proportion of coding neuron')
        text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    end
    min_max_list=[min_max_list;2,ii,r,p];
end
%%
if isequal(feat_reg_map,sens_reg_map)
    s_list=sortrows(min_max_list,3,'descend');
    s_list=s_list(isfinite(s_list(:,3)),:);
    fhb=figure('Color','w','Position',[32,32,1200,225]);
    hold on
    bhto=bar(find(s_list(:,1)==1),s_list(s_list(:,1)==1,3),0.6,"grouped",'white');
    bhfrom=bar(find(s_list(:,1)==2),s_list(s_list(:,1)==2,3),0.6,"grouped",'black');
    ylabel('Connectivity-coding proportion correlation (Pearson''s r)')
    xlbl=cell(size(s_list,1),1);
    xlbl(s_list(:,1)==1)=idmap.ccfid2reg.values(num2cell(sink_ccfid(s_list(s_list(:,1)==1,2))));
    xlbl(s_list(:,1)==2)=idmap.ccfid2reg.values(num2cell(src_ccfid(s_list(s_list(:,1)==2,2))));
    xlbl=cellfun(@(x) x{1}, xlbl,'UniformOutput',false);
    set(gca(),'XTick',1:size(s_list,1),'XTickLabel',xlbl);
    legend([bhto,bhfrom],{'Connectivity to this region','Connectivity from this region'})
    exportgraphics(fhb,'Sense_proportion_connectivity_corr_bars.pdf','ContentType','vector');
    [maxr,maxidx]=max(min_max_list(:,3));
    [minr,minidx]=min(min_max_list(:,3));

    % hand picked as guided by max-min data
    idmap.ccfid2reg(sink_ccfid(20)) %AON
    idmap.ccfid2reg(sink_ccfid(50)) %GPe

    allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
    intersect_regs=intersect(allen_regs,chregs);
    idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
    feat_prop_cell=feat_reg_map.values(intersect_regs);
    feat_prop=cellfun(@(x) x(1),feat_prop_cell);

    AON_density=log10(src_sink_mat(20,idx4corr)); %AON
    GPe_density=log10(src_sink_mat(50,idx4corr)); %GPe
    allen_density=AON_density-GPe_density;
    [r,p]=corr(allen_density.',feat_prop);

    fhs=figure('Color','w');
    scatter(allen_density,feat_prop,4,'red','filled','o')
    text(allen_density,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
    xlabel('(To-AON) - (To-GPe) Projection density (log10)')
    ylabel('Proportion of coding neuron')
    title(sprintf('r=%.3f, p=%.3f',r,p));
    exportgraphics(fhs,'Sense_proportion_connectivity_corr_scatter.pdf','ContentType','vector');
end

%%
if isequal(feat_reg_map,dur_reg_map)
    [maxr,maxidx]=max(min_max_list(:,3));
    [minr,minidx]=min(min_max_list(:,3));

    % hand picked as guided by max-min data
    idmap.ccfid2reg(sink_ccfid(45)) %AAA
    idmap.ccfid2reg(src_ccfid(8)) %ACA

    allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
    intersect_regs=intersect(allen_regs,chregs);
    idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
    feat_prop_cell=feat_reg_map.values(intersect_regs);
    feat_prop=cellfun(@(x) x(1),feat_prop_cell);

    AAA_density=log10(src_sink_mat(45,idx4corr)); %AAA
    ACA_density=log10(src_sink_mat(idx4corr,8)); %ACA
    allen_density=AON_density-GPe_density;
    [r,p]=corr(allen_density.',feat_prop);

    fhs=figure('Color','w');
    scatter(allen_density,feat_prop,4,'red','filled','o')
    text(allen_density,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
    xlabel('AAA-ACA Projection density (log10)')
    ylabel('Proportion of coding neuron')
    text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    
end
