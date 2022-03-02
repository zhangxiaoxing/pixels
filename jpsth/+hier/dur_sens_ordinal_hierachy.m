%% CONST

dur_ord_bin=1:17%[1:5,10:17];

meta=ephys.util.load_meta();
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
% ctxsel=strcmp(meta.reg_tree(2,:),'CTX') & ~strcmp(meta.reg_tree(5,:),'');
% cnusel=strcmp(meta.reg_tree(2,:),'CNU') & ~strcmp(meta.reg_tree(5,:),'');
% chregs=unique(meta.reg_tree(5,ctxsel|cnusel));
BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));

cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
grey_regs=grey_regs(cnt>100);

% preferably load disk data

%% Sense|odor
waveid=ephys.get_wave_id(meta.sess,meta.allcid);
sens_reg_map=containers.Map();
both_reg_map=containers.Map();
for r=grey_regs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    pos=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & waveid>0);
    both_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & waveid>4);
    sens_reg_map(r{1})=[pos/cnt,pos,cnt];
    both_reg_map(r{1})=[both_cnt/cnt,pos,cnt];
end

%% duration
dur=wave.duration_wave('quick_merge',true,'align_test',false,'ctx',false);
dur_reg_map=containers.Map();
for r=grey_regs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    pos=nnz(strcmp(meta.reg_tree(5,:),r{1}) & any((dur.wrsp(:,dur_ord_bin)<0.05).'));
    dur_reg_map(r{1})=[pos/cnt,pos,cnt];
end

%% ordinal
ord=load('Ordinal_MI_Combined_WelltrainedTrials_Corr.mat','out_');
ord_reg_map=containers.Map();
avai_sess=unique(ord.out_.su_meta(:,1));
sess_sel=ismember(meta.sess,avai_sess);
for r=grey_regs
    cnt=nnz(strcmp(meta.reg_tree(5,sess_sel),r{1}));
    pos=nnz(strcmp(meta.reg_tree(5,sess_sel),r{1}) & any((ord.out_.anovanp(:,dur_ord_bin)<0.05).'));
    ord_reg_map(r{1})=[pos/cnt,pos,cnt];
end

%% assert dimension equality
if false % one-time consistency check
    assert(isequal(size(dur.wrsp),[33028,17]))

    assert(isequal(meta.sess,dur.meta(:,1)));
    assert(isequal(meta.allcid,dur.meta(:,2)));

    assert(isequal(meta.sess(sess_sel),ord.out_.su_meta(:,1)));
    assert(isequal(meta.allcid(sess_sel),ord.out_.su_meta(:,3)));
end

%% allen data
sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');

nan_rmv=any(~isfinite(sink_src_mat),2);
sink_ccfid=sink_ccfid(~nan_rmv);
sink_src_mat=sink_src_mat(~nan_rmv,:);

src_idx_map=containers.Map(src_ccfid,1:numel(src_ccfid)); %{ccfid:proj_dense_mat_idx}
sink_idx_map=containers.Map(sink_ccfid,1:numel(sink_ccfid)); %{ccfid:proj_dense_mat_idx}

if false
    figure();
    imagesc(sink_src_mat);
end

%%
out=struct();
for type="both"%["sense","both","dur","ord"]
    switch type
        case 'sense'
            feat_reg_map=sens_reg_map;
        case 'both'
            feat_reg_map=both_reg_map;            
        case 'dur'
            feat_reg_map=dur_reg_map;
        case 'ord'
            feat_reg_map=ord_reg_map;
    end

    %% To GLM

    if false %to_GLM
        wave.connectivity_proportion_GLM
    else
        min_max_list=[];
        %% efferent_proj_dense(src) ~ feature_proportion
        allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
        intersect_regs=intersect(allen_regs,grey_regs);
        idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
        feat_prop_cell=feat_reg_map.values(intersect_regs);
        feat_prop=cellfun(@(x) x(1),feat_prop_cell);

        for ii=1:numel(sink_ccfid)
            allen_density=log10(sink_src_mat(ii,idx4corr)); % from idx4corr, to one alternating target
            [r,p]=corr(allen_density.',feat_prop);
            %         if p<0.05 && false
            %             figure('Color','w')
            %             scatter(allen_density,feat_prop,4,'red','filled','o')
            %             text(allen_density,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
            %             title(['Allen connectivity from each region to', idmap.ccfid2full(sink_ccfid(ii))])
            %             xlabel('Projection density (log10)')
            %             ylabel('Proportion of coding neuron')
            %             text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
            %         end
            min_max_list=[min_max_list;1,ii,r,p];
        end

        %% afferent_proj_dense(sink) ~ feature_proportion

        allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(sink_ccfid)),'UniformOutput',false);
        intersect_regs=intersect(allen_regs,grey_regs);
        idx4corr=cell2mat(sink_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
        feat_prop_cell=feat_reg_map.values(intersect_regs);
        feat_prop=cellfun(@(x) x(1),feat_prop_cell);

        for ii=1:numel(src_ccfid)
            allen_density=log10(sink_src_mat(idx4corr,ii)); % from one alternating target, to idx4corr
            [r,p]=corr(allen_density,feat_prop);
            %         if p<0.05 && false
            %             figure('Color','w')
            %             scatter(allen_density,feat_prop,4,'red','filled','o')
            %             text(allen_density,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
            %             title(['Allen connectivity to each region from', idmap.ccfid2full(src_ccfid(ii))])
            %             xlabel('Projection density (log10)')
            %             ylabel('Proportion of coding neuron')
            %             text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
            %         end
            min_max_list=[min_max_list;2,ii,r,p];
        end
        %%
        s_list=sortrows(min_max_list,3,'descend');
        s_list=s_list(isfinite(s_list(:,3)),:);

        xlbl=cell(size(s_list,1),1);
        xlbl(s_list(:,1)==1)=idmap.ccfid2reg.values(num2cell(sink_ccfid(s_list(s_list(:,1)==1,2))));
        xlbl(s_list(:,1)==2)=idmap.ccfid2reg.values(num2cell(src_ccfid(s_list(s_list(:,1)==2,2))));
        xlbl=cellfun(@(x) x{1}, xlbl,'UniformOutput',false);

        out.(type).r_p=s_list;
        out.(type).reg=xlbl;
        fhb=figure('Color','w','Position',[32,32,1200,225]);
        hold on
        bhto=bar(find(s_list(:,1)==1),s_list(s_list(:,1)==1,3),0.6,"grouped",'white');
        bhfrom=bar(find(s_list(:,1)==2),s_list(s_list(:,1)==2,3),0.6,"grouped",'black');
        ylabel('Connectivity-coding proportion correlation (Pearson''s r)')
        set(gca(),'XTick',1:size(s_list,1),'XTickLabel',xlbl);
        legend([bhto,bhfrom],{'Connectivity to this region','Connectivity from this region'})

        if false
            if isequal(feat_reg_map,sens_reg_map)
                exportgraphics(fhb,'Sense_proportion_connectivity_corr_bars.pdf','ContentType','vector');
                %for next graph
                idmap.ccfid2reg(sink_ccfid(20)) %AON
                idmap.ccfid2reg(sink_ccfid(50)) %GPe

                POS_density=log10(sink_src_mat(20,idx4corr)); %AON
                Neg_density=log10(sink_src_mat(50,idx4corr)); %GPe

            elseif isequal(feat_reg_map,dur_reg_map)
                exportgraphics(fhb,'Dur_proportion_connectivity_corr_bars.pdf','ContentType','vector');
                %for next graph
                idmap.ccfid2reg(src_ccfid(17)) %AON
                idmap.ccfid2reg(src_ccfid(6)) % AUD

                POS_density=log10(sink_src_mat(idx4corr,17)); %AON
                Neg_density=log10(sink_src_mat(idx4corr,6)); %AUD
            end

            % hand picked as guided by max-min data
            allen_density=reshape(POS_density-Neg_density,[],1);
            [r,p]=corr(allen_density,feat_prop);

            fhs=figure('Color','w');
            scatter(allen_density,feat_prop,4,'red','filled','o')
            text(allen_density,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
            ylabel('Proportion of coding neuron')
            title(sprintf('r=%.3f, p=%.3f',r,p));
            if isequal(feat_reg_map,sens_reg_map)
                xlabel('(To-AON) - (To-GPe) Projection density (log10)')
                exportgraphics(fhs,'Sense_proportion_connectivity_corr_scatter.pdf','ContentType','vector');
            elseif isequal(feat_reg_map,dur_reg_map)
                xlabel('(From-AON) - (From-AUD) Projection density (log10)')
                exportgraphics(fhs,'Dur_proportion_connectivity_corr_scatter.pdf','ContentType','vector');
            end

            %%
            if isequal(feat_reg_map,dur_reg_map)
                [maxr,maxidx]=max(min_max_list(:,3));
                [minr,minidx]=min(min_max_list(:,3));

                % hand picked as guided by max-min data
                idmap.ccfid2reg(sink_ccfid(45)) %AAA
                idmap.ccfid2reg(src_ccfid(8)) %ACA

                allen_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
                intersect_regs=intersect(allen_regs,grey_regs);
                idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
                feat_prop_cell=feat_reg_map.values(intersect_regs);
                feat_prop=cellfun(@(x) x(1),feat_prop_cell);

                AAA_density=log10(sink_src_mat(45,idx4corr)); %AAA
                ACA_density=log10(sink_src_mat(idx4corr,8)); %ACA
                allen_density=POS_density-Neg_density;
                [r,p]=corr(allen_density.',feat_prop);

                fhs=figure('Color','w');
                scatter(allen_density,feat_prop,4,'red','filled','o')
                text(allen_density,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
                xlabel('AAA-ACA Projection density (log10)')
                ylabel('Proportion of coding neuron')
                text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');

            end
        end
    end
end
% 3 feature 3d-PCA
if false
    sinkregs=intersect(intersect(out.sense.reg(out.sense.r_p(:,1)==1),out.dur.reg(out.dur.r_p(:,1)==1)),out.ord.reg(out.ord.r_p(:,1)==1));
    srcregs=intersect(intersect(out.sense.reg(out.sense.r_p(:,1)==2),out.dur.reg(out.dur.r_p(:,1)==2)),out.ord.reg(out.ord.r_p(:,1)==2));

    src_sink_regs=intersect(sinkregs,srcregs);

    fh=figure('Color','w');
    hold on;

    for r=reshape(sinkregs,1,[])
        rr=r{1};
        rh=plot3(out.sense.r_p(strcmp(out.sense.reg,rr) & out.sense.r_p(:,1)==1,3),...
            out.dur.r_p(strcmp(out.dur.reg,rr) & out.dur.r_p(:,1)==1,3),...
            out.ord.r_p(strcmp(out.ord.reg,rr) & out.ord.r_p(:,1)==1,3),'ro');
    end

    for r=reshape(srcregs,1,[])
        rr=r{1};
        bh=plot3(out.sense.r_p(strcmp(out.sense.reg,rr) & out.sense.r_p(:,1)==2,3),...
            out.dur.r_p(strcmp(out.dur.reg,rr) & out.dur.r_p(:,1)==2,3),...
            out.ord.r_p(strcmp(out.ord.reg,rr) & out.ord.r_p(:,1)==2,3),'bo')
    end

    for r=reshape(src_sink_regs,1,[])
        rr=r{1};
        plot3([out.sense.r_p(strcmp(out.sense.reg,rr) & out.sense.r_p(:,1)==1,3),out.sense.r_p(strcmp(out.sense.reg,rr) & out.sense.r_p(:,1)==2,3)],...
            [out.dur.r_p(strcmp(out.dur.reg,rr) & out.dur.r_p(:,1)==1,3),out.dur.r_p(strcmp(out.dur.reg,rr) & out.dur.r_p(:,1)==2,3)],...
            [out.ord.r_p(strcmp(out.ord.reg,rr) & out.ord.r_p(:,1)==1,3),out.ord.r_p(strcmp(out.ord.reg,rr) & out.ord.r_p(:,1)==2,3)],'-','Color',[0.5,0.5,0.5])
    end
    v=VideoWriter('Sens_dur_ord_connectivity_proportion_corr.mp4');
    open(v)
    for az=-130:360-130
        view([az,45])
        writeVideo(v,getframe(gcf()))
    end
    close(v)
end
