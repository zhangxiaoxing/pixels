
%% CONST
corr1=true;
plot1=true;
corr2=false;
plot2=false;
if isunix
    load('map_cells.mat','map_cells');
    maxNumCompThreads(64)
else
    load('map_cells.mat','map_cells');
    maxNumCompThreads(4)
end

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
grey_regs=ephys.getGreyRegs();

sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');
% grey_regs assume in work space
src_idx_map=containers.Map(src_ccfid,1:numel(src_ccfid)); %{ccfid:proj_dense_mat_idx}
sink_idx_map=containers.Map(sink_ccfid,1:numel(sink_ccfid)); %{ccfid:proj_dense_mat_idx}

allen_src_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
intersect_regs=intersect(allen_src_regs,grey_regs);


%% -> feature_region_map entry point
one_reg_corr_list=[];
two_reg_corr_list=[];
for featii=1:4 %mixed sample, exclusive sample, mixed duration, exclusive duration
    for epochii=1:3 %sample,early delay, late delay, test, post_reward, pre_next_sample

        idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
        feat_reg_map=map_cells{epochii,featii};

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
        if corr1
            for regii=1:size(glmxmat,1)
                [r,p]=corr(glmxmat(regii,:).',feat_prop);
                one_reg_corr_list=[one_reg_corr_list;featii,epochii,regii,double(glmxmeta(regii,:)),r,p];
                %====================================^^^1^^^^^^2^^^^^^3^^^^^^^^^^^4,5,6^^^^^^^^^^^^^7^8^^
            end
            save('one_reg_corr_list.mat','one_reg_corr_list')

            %
            if plot1 % plot one-region: coding proportion correlation bars
                s_list=sortrows(one_reg_corr_list(one_reg_corr_list(:,4)==2 & one_reg_corr_list(:,1)==featii & one_reg_corr_list(:,2)==epochii,:),7,'descend');
%                 s_list=s_list(isfinite(s_list(:,7)),:);
                s_list=s_list(1:15,:);
                xlbl=cell(size(s_list,1),1);
                xlbl(s_list(:,4)==1)=idmap.ccfid2reg.values(num2cell(s_list(s_list(:,4)==1,6)));
                xlbl(s_list(:,4)==2)=idmap.ccfid2reg.values(num2cell(s_list(s_list(:,4)==2,6)));
                xlbl=cellfun(@(x) x{1}, xlbl,'UniformOutput',false);

                fhb=figure('Color','w','Position',[32,32,215,225]);
                hold on
%                 bhto=bar(find(s_list(:,4)==1),s_list(s_list(:,4)==1,7),0.6,'white');
%                 bhfrom=bar(find(s_list(:,4)==2),s_list(s_list(:,4)==2,7),0.6,'black');
                bh=bar(s_list(:,7),0.6,'FaceColor','flat');
                bh.CData(s_list(:,4)==1,:)=repmat([1,1,1],nnz(s_list(:,4)==1),1);
                bh.CData(s_list(:,4)==2,:)=repmat([0,0,0],nnz(s_list(:,4)==2),1);
                ylabel('Connectivity-coding proportion correlation (Pearson''s r)')
                set(gca(),'XTick',1:size(s_list,1),'XTickLabel',xlbl,'XTickLabelRotation',90);
                title(sprintf('epoch%d-feat%d',epochii,featii))
%                 legend([bhto,bhfrom],{'Connectivity to this region','Connectivity from this region'})
            end
        end
        %% two region glm interaction
        if corr2
            comb2=nchoosek(1:size(glmxmat,1),2);
            mdlid_rsq_AICC=[];
            for jj=1:size(comb2,1)
                if rem(jj,1000)==0
                    disp(jj)
                end
                allen_A=glmxmat(comb2(jj,1),:).'; % from one alternating target, to idx4corr
                allen_B=glmxmat(comb2(jj,2),:).';
                mdl=fitglm([allen_A,allen_B],feat_prop,'interactions'); % (1|2)jj => glmxmat==glmxmeta => (sink_ccfid|src_ccfid)
                two_reg_corr_list=[two_reg_corr_list;featii,epochii,jj,double(glmxmeta(comb2(jj,1),:)),double(glmxmeta(comb2(jj,2),:)),mdl.Coefficients.Estimate.',mdl.Rsquared.Ordinary,mdl.ModelCriterion.AICc];
                %============================================^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                mdlid_rsq_AICC=[mdlid_rsq_AICC;jj,mdl.Rsquared.Ordinary,mdl.ModelCriterion.AICc];
                %     if mdl.Rsquared.Ordinary>maxrsq
                %         maxrsq=mdl.Rsquared.Ordinary;
                %         maxidx=jj;
                %     end
            end
            if plot2
                mdlid_rsq_AICC=sortrows(mdlid_rsq_AICC,3);
                % keyboard();
                ytk=cell(0,0);
                for kk=1:20
                    maxidx=mdlid_rsq_AICC(kk,1);
                    if glmxmeta(comb2(maxidx,1),1)==1
                        reg1=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,1),2)));
                        reg1=['To ',reg1{1}];
                    else
                        reg1=(idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(maxidx,1),2))));
                        reg1=['From ',reg1{1}];
                    end
                    if glmxmeta(comb2(maxidx,2),1)==1
                        reg2=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,2),2)));
                        reg2=['To ',reg2{1}];
                    else
                        reg2=idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(maxidx,2),2)));
                        reg2=['From ',reg2{1}];
                    end
                    disp([reg1,'-',reg2,' ',num2str(mdlid_rsq_AICC(kk,2:3))])
                    ytk{end+1}=[reg1,'-',reg2];
                end

                fh=figure('Color','w','Position',[32,32,400,800]);
                bar(sqrt(mdlid_rsq_AICC(1:10,2)),'Horizontal','on')
                set(gca(),'YDir','reverse','YTick',1:10,'YTickLabel',ytk(1:10))
                xlabel('Person''s r')
                xlim([0,1])
                title(sprintf('selectivity %d - epoch %d',featii,epochii));
                exportgraphics(fh,sprintf('frac_allen_mdls_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')
                close(fh)

                maxidx=mdlid_rsq_AICC(1,1);
                allen_A=glmxmat(comb2(maxidx,1),:).'; % from one alternating target, to idx4corr
                allen_B=glmxmat(comb2(maxidx,2),:).';
                mdl=fitglm([allen_A,allen_B],feat_prop,'interactions');


                disp(sqrt(mdl.Rsquared.Ordinary));
                if glmxmeta(comb2(maxidx,1),1)==1
                    reg1=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,1),2)));
                    reg1=['To_',reg1{1}];
                else
                    reg1=(idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(maxidx,1),2))));
                    reg1=['From_',reg1{1}];
                end
                if glmxmeta(comb2(maxidx,2),1)==1
                    reg2=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,2),2)));
                    reg2=['To_',reg2{1}];
                else
                    reg2=idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(maxidx,2),2)));
                    reg2=['From_',reg2{1}];
                end

                fh=figure('Color','w','Position',[32,32,320,320]);
                hold on
                for ll=1:numel(feat_prop)
                    c=ephys.getRegColor(intersect_regs{ll},'large_area',true);
                    scatter(mdl.Fitted.Response(ll),feat_prop(ll),4,c,'filled','o')
                    text(mdl.Fitted.Response(ll),feat_prop(ll),intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
                end
                title(sprintf('Coding proportion ~ %.3f [%s] + %.3f [%s] +%.3f [%s]', ...
                    mdl.Coefficients.Estimate(2),reg1,mdl.Coefficients.Estimate(3),reg2,mdl.Coefficients.Estimate(4),'interaction'), ...
                    'Interpreter','none')
                xlabel(sprintf('Model %d prediction',maxidx))
                ylabel('Proportion of coding neuron')
                set(gca,'XScale','log','YScale','log')
                text(min(xlim()),max(ylim()),sprintf(' r = %.3f, AIC = %.1f',sqrt(mdl.Rsquared.Ordinary),mdl.ModelCriterion.AIC),'HorizontalAlignment','left','VerticalAlignment','top');
                % xlim([0.15,0.5])
                % ylim([0.15,0.5])
                exportgraphics(fh,sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')
                close(fh)
                %         keyboard()

            end
            save('two_reg_corr_list.mat','two_reg_corr_list')
        end
    end
end

return
%%
comb3=nchoosek(1:size(glmxmat,1),3);
max3rsqp=0;
max3idxp=-1;
for jj=127706%1:size(comb3,1)
    if rem(jj,1000)==0
        disp(jj)
    end
    allen_A=glmxmat(comb3(jj,1),:).'; % from one alternating target, to idx4corr
    allen_B=glmxmat(comb3(jj,2),:).';
    allen_C=glmxmat(comb3(jj,3),:).';
    mdl=fitglm([allen_A,allen_B,allen_C],feat_prop,'poly111');
    if mdl.Rsquared.Ordinary>max3rsqp
        max3rsqp=mdl.Rsquared.Ordinary;
        max3idxp=jj;
    end

    if mdl.Rsquared.Ordinary>0.64
        disp(jj)
        disp(sqrt(mdl.Rsquared.Ordinary));
        if glmxmeta(comb3(jj,1),1)==1
            reg1=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb3(jj,1),2)));
            reg1=['To_',reg1{1}];
        else
            reg1=(idmap.ccfid2reg(src_ccfid(glmxmeta(comb3(jj,1),2))));
            reg1=['From_',reg1{1}];
        end
        if glmxmeta(comb3(jj,2),1)==1
            reg2=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb3(jj,2),2)));
            reg2=['To_',reg2{1}];
        else
            reg2=idmap.ccfid2reg(src_ccfid(glmxmeta(comb3(jj,2),2)));
            reg2=['From_',reg2{1}];
        end

        if glmxmeta(comb3(jj,3),1)==1
            reg3=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb3(jj,3),2)));
            reg3=['To_',reg3{1}];
        else
            reg3=idmap.ccfid2reg(src_ccfid(glmxmeta(comb3(jj,3),2)));
            reg3=['From_',reg3{1}];
        end

        fh=figure('Color','w','Position',[32,32,900,650]);
        scatter(mdl.Fitted.Response,feat_prop,4,'red','filled','o')
        text(mdl.Fitted.Response,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
        title(sprintf('Coding proportion ~ %.3f [%s] + %.3f [%s] +%.3f [%s]', ...
            mdl.Coefficients.Estimate(2),reg1,mdl.Coefficients.Estimate(3),reg2,mdl.Coefficients.Estimate(4),reg3), ...
            'Interpreter','none')
        xlabel(sprintf('Model %d prediction',jj))
        ylabel('Proportion of coding neuron')
        text(min(xlim()),max(ylim()),sprintf(' r = %.3f, AIC = %.1f',sqrt(mdl.Rsquared.Ordinary),mdl.ModelCriterion.AIC),'HorizontalAlignment','left','VerticalAlignment','top');
        %         keyboard()
    end
end

