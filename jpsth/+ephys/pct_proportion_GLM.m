% Likely obsolete as of 2023.03.08
% pending delete

function fh=pct_proportion_GLM(map_cells,corr_type,opt)
arguments
    map_cells (1,:) cell
    corr_type (1,:) char {mustBeMember(corr_type,{'Spearman','PearsonLogLog','PearsonLinearLog'})}
%     opt.skip_efferent (1,1) logical = true
%     opt.corr1 (1,1) logical =true;
%     opt.plot1 (1,1) logical =true;
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    opt.stats_type (1,:) char
    opt.data_type (1,:) char
    opt.feat_tag
    
end

if ~isfield(opt,'feat_tag') || isempty(opt.feat_tag)
    warning('Missing data tag');
end

%map_cells from K:\code\jpsth\+ephys\Both_either_reg_bars.m
%% allen ref value
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
grey_regs=ephys.getGreyRegs('range',opt.range);

sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');
src_idx_map=containers.Map(src_ccfid,1:numel(src_ccfid)); %{ccfid:proj_dense_mat_idx}
sink_idx_map=containers.Map(sink_ccfid,1:numel(sink_ccfid)); %{ccfid:proj_dense_mat_idx}

allen_src_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
map_regs=intersect(allen_src_regs,grey_regs);

% disp({opt.range,numel(map_regs)})

%% -> feature_region_map entry point
for featii=1:numel(map_cells) % both, either, summed
    one_reg_corr_list=[];
    feat_reg_map=map_cells{featii};
    intersect_regs=intersect(map_regs,feat_reg_map.keys());
    idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));

    feat_prop_cell=feat_reg_map.values(intersect_regs);
    feat_prop=cellfun(@(x) x(1),feat_prop_cell);

    [glmxmat,glmxmeta]=deal([]);
    % efferent_proj_dense(soruce) ~ feature_proportion
    for ii=1:numel(sink_ccfid)
        %     allen_density=log10(sink_src_mat(ii,idx4corr)+1e-12); % from idx4corr, to one alternating target
        allen_density=sink_src_mat(ii,idx4corr);
        if nnz(allen_density~=0)<10
            continue
        end
        glmxmat=[glmxmat;allen_density];
        glmxmeta=[glmxmeta;1,ii,sink_ccfid(ii)];
    end

    % afferent_proj_dense(sink) ~ feature_proportion
    idx4corr=cell2mat(sink_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
    for ii=1:numel(src_ccfid)
        allen_density=sink_src_mat(idx4corr,ii); % from one alternating target, to idx4corr
        if nnz(allen_density~=0)<10
            disp(src_ccfid)
            continue
        end
        glmxmat=[glmxmat;allen_density.'];
        glmxmeta=[glmxmeta;2,ii,src_ccfid(ii)];
    end

    %% one region corr, bar plot
%     if opt.corr1
        for regii=1:size(glmxmat,1)
            switch(corr_type)
                case 'Spearman'
                    [r,p]=corr(glmxmat(regii,:).',feat_prop,'type','Spearman');
                case 'PearsonLogLog'
                    vsel=(glmxmat(regii,:).')>0 & feat_prop>0;
                    [r,p]=corr(reshape(log10(glmxmat(regii,vsel)),[],1),...
                        reshape(log10(feat_prop(vsel)),[],1),'type','Pearson');
                case 'PearsonLinearLog'
                    vsel=(glmxmat(regii,:).')>0;
                    [r,p]=corr(reshape(log10(glmxmat(regii,vsel)),[],1),...
                        reshape(feat_prop(vsel),[],1),'type','Pearson');
            end
            
            
            one_reg_corr_list=[one_reg_corr_list;featii,NaN,regii,double(glmxmeta(regii,:)),r,p];
            %====================================^^^1^^^^2^^^^3^^^^^^^^^^^4,5,6^^^^^^^^^^^^^7^8^^
        end
%         save('one_reg_corr_list.mat','one_reg_corr_list')

        %
%         if opt.plot1 % plot one-region: coding proportion correlation bars
            s_list=sortrows(one_reg_corr_list(one_reg_corr_list(:,1)==featii,:),7,'descend');
            %                 s_list=s_list(isfinite(s_list(:,7)),:);
%             if opt.skip_efferent
                s_list=s_list(s_list(:,4)==2,:);
%             end
            
            s_list=s_list([1:5,end-4:end],:);
            xlbl=idmap.ccfid2reg.values(num2cell(s_list(:,6)));
            xlbl=cellfun(@(x) x{1}, xlbl,'UniformOutput',false);
            
            % model list
            fh.(sprintf('feat%d',featii))=figure('Color','w','Position',[32,32,650,700]);
            tiledlayout(2,5);
            nexttile(4,[1,2])
            hold on
            %                 bhto=bar(find(s_list(:,4)==1),s_list(s_list(:,4)==1,7),0.6,'white');
            %                 bhfrom=bar(find(s_list(:,4)==2),s_list(s_list(:,4)==2,7),0.6,'black');
            bh=bar(s_list(:,7),0.8,'FaceColor','flat','Horizontal','on');
            bh.CData(s_list(:,4)==1,:)=repmat([0.5,0.5,0.5],nnz(s_list(:,4)==1),1);
            bh.CData(s_list(:,4)==2,:)=repmat([0,0,0],nnz(s_list(:,4)==2),1);
            xlabel('Connectivity-coding proportion correlation (Pearson''s r)')
            set(gca(),'YTick',1:size(s_list,1),'YTickLabel',xlbl,'YDir','reverse');
            if isfield(opt,'feat_tag') && ~isempty(opt.feat_tag)
                title(opt.feat_tag(featii));
            else
                title(sprintf('feature%d',featii))
            end
            xlim([-1,1]);
            %                 legend([bhto,bhfrom],{'Connectivity to this region','Connectivity from this region'})
%             exportgraphics(fh.bar1,sprintf('frac_allen_mdls_selec%d.pdf',featii),'ContentType','vector')
            % positive top
            nexttile(1,[1,3])
            hold on
            glmidx=find(glmxmeta(:,1)==s_list(1,4) & glmxmeta(:,2)==s_list(1,5));
            for ll=1:numel(feat_prop)
                c=ephys.getRegColor(intersect_regs{ll},'large_area',true);
                scatter(glmxmat(glmidx,ll),feat_prop(ll),4,c,'filled','o')
                text(glmxmat(glmidx,ll),feat_prop(ll),intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
            end
            xlabel('Projection density (px/px)');
            ylabel('Regional averaged feature index')
            if strcmp(corr_type,'PearsonLinearLog')
                set(gca,'XScale','log','YScale','linear')
            else
                set(gca,'XScale','log','YScale','log')
            end
            title(sprintf(' r = %.3f, p = %.3f',s_list(1,7),s_list(1,8)));

            nexttile(6,[1,3])
            hold on
            glmidx=find(glmxmeta(:,1)==s_list(end,4) & glmxmeta(:,2)==s_list(end,5));
            for ll=1:numel(feat_prop)
                c=ephys.getRegColor(intersect_regs{ll},'large_area',true);
                scatter(glmxmat(glmidx,ll),feat_prop(ll),4,c,'filled','o')
                text(glmxmat(glmidx,ll),feat_prop(ll),intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
            end
            xlabel('Projection density (px/px)');
            ylabel('Regional averaged feature index')
            if strcmp(corr_type,'PearsonLinearLog')
                set(gca,'XScale','log','YScale','linear')
            else
                set(gca,'XScale','log','YScale','log')
            end
            title(sprintf(' r = %.3f, p = %.3f',s_list(end,7),s_list(end,8)));

            th=nexttile(9,[1,2]);
            tbl=cell(0);
            if isfield(opt,'stats_type') && ~isempty(opt.stats_type)
                tbl=[tbl;'Stats';opt.stats_type];
            end
            if isfield(opt,'data_type') && ~isempty(opt.data_type)
                tbl=[tbl;'Data';opt.data_type];
            end
            tbl=[tbl;'Range';opt.range];
            ephys.util.figtable(gcf(),th,tbl)

%             exportgraphics(fh.corr1,sprintf('frac_allen_scatter_selec%d.pdf',featii),'ContentType','vector')

 %         end
%     end
end