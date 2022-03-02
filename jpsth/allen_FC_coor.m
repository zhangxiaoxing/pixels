idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
[sig,pair]=bz.load_sig_pair('pair',true,'CTX',false);


sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/CH_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/CH_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');

nan_rmv=any(~isfinite(sink_src_mat),2);
sink_ccfid=sink_ccfid(~nan_rmv);
sink_src_mat=sink_src_mat(~nan_rmv,:);

src_idx_map=containers.Map(src_ccfid,1:numel(src_ccfid)); %{ccfid:proj_dense_mat_idx}
sink_idx_map=containers.Map(sink_ccfid,1:numel(sink_ccfid)); %{ccfid:proj_dense_mat_idx}

fc_reg_ids=unique(pair.reg(:,5,1));

corrs=[];
skipcnt=0;
for srcid=reshape(src_ccfid,1,[])
    for sinkid=reshape(sink_ccfid,1,[])
        allenPD=sink_src_mat(sink_ccfid==sinkid,src_ccfid==srcid);
        paircnt=nnz(pair.reg(:,5,1)==srcid & pair.reg(:,5,2)==sinkid);
        if paircnt==0
            disp([idmap.ccfid2full(srcid),idmap.ccfid2full(sinkid)])
            skipcnt=skipcnt+1;
        end
        sigcnt=nnz(sig.reg(:,5,1)==srcid & sig.reg(:,5,2)==sinkid);
        corrs=[corrs;double(srcid),double(sinkid),allenPD,log10(allenPD),sigcnt,paircnt,sigcnt/paircnt];
    end
end

corrs_raw=corrs;
corrs=corrs_raw(corrs_raw(:,6)>=1000,:);
corrs(:,7)=corrs(:,7)+1e-12;
unnest=@(x) cellfun(@(y) y{1}, x,'UniformOutput',false);

unnest6=@(x) cellfun(@(y) y{6}, x,'UniformOutput',false);

src_dep5=unnest6(idmap.reg2tree.values((unnest(idmap.ccfid2reg.values(num2cell(corrs(:,1)))))));
sink_dep5=unnest6(idmap.reg2tree.values((unnest(idmap.ccfid2reg.values(num2cell(corrs(:,2)))))));

% figure()
% scatter(corrs(:,3),corrs(:,7))
% 
% figure()
% scatter(corrs(:,4),corrs(:,7))

srcsel=strcmp(src_dep5,'OLF') & ~strcmp(sink_dep5,'OLF');
sinksel=~strcmp(src_dep5,'OLF') & strcmp(sink_dep5,'OLF');
bothsel=strcmp(src_dep5,'OLF') & strcmp(sink_dep5,'OLF');
othersel=~(srcsel|sinksel|bothsel);
figure()
hold on
scatter(corrs(srcsel,4),log10(corrs(srcsel,7)),9,'r','filled')
scatter(corrs(sinksel,4),log10(corrs(sinksel,7)),9,'b','filled')
scatter(corrs(bothsel,4),log10(corrs(bothsel,7)),9,'m','filled')
scatter(corrs(othersel,4),log10(corrs(othersel,7)),9,'k','filled')

[r,p]=corr(corrs(:,3),corrs(:,7))

% [r,p]=corr(corrs(:,3),corrs(:,7))