function [fh,hiermap]=connectivity_proportion_GLM(map_cells,corr_type,opt)
arguments
    map_cells (1,:) struct
    corr_type (1,:) char {mustBeMember(corr_type,{'Spearman','PearsonLogLog','PearsonLinearLog'})}
    opt.skip_efferent (1,1) logical = true
    opt.corr1 (1,1) logical =true;
    opt.plot1 (1,1) logical =true;
    opt.corr2 (1,1) logical =false;
    opt.plot2 (1,1) logical =false;
    opt.corr3 (1,1) logical =false;
    opt.plot3 (1,1) logical =false;
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    opt.stats_type (1,:) char
    opt.data_type (1,:) char
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end

if opt.corr2 || opt.corr3
    error("Need to remove injection site from correlation, not done yet")
end


%map_cells from K:\code\jpsth\+ephys\Both_either_reg_bars.m

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
grey_regs=ephys.getGreyRegs('range',opt.range,'criteria',opt.criteria);

sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');
src_idx_map=containers.Map(src_ccfid,1:numel(src_ccfid)); %{ccfid:proj_dense_mat_idx}
sink_idx_map=containers.Map(sink_ccfid,1:numel(sink_ccfid)); %{ccfid:proj_dense_mat_idx}

allen_src_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
map_regs=intersect(allen_src_regs,grey_regs);

% disp({opt.range,numel(map_regs)})


%% -> feature_region_map entry point
epochii=1;
fns=fieldnames(map_cells);
for featii=1:numel(fns) % both, either, summed
    one_reg_corr_list=[];
%     two_reg_corr_list=[];
    feat_reg_map=map_cells.(fns{featii});
    intersect_regs=intersect(map_regs,feat_reg_map.keys());
    idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));

    feat_prop=cell2mat(feat_reg_map.values(intersect_regs));
    if size(feat_prop,2)>1
        feat_prop=feat_prop(:,1);
    end
    [allen_mat,allen_meta]=deal([]);
    % efferent_proj_dense(soruce) ~ feature_proportion
    if ~opt.skip_efferent
        for ii=1:numel(sink_ccfid)
            %     allen_density=log10(sink_src_mat(ii,idx4corr)+1e-12); % from idx4corr, to one alternating target
            allen_density=sink_src_mat(ii,idx4corr);
            if nnz(allen_density~=0)<10
                continue
            elseif ~ismember(idmap.ccfid2reg(sink_ccfid(ii)),map_regs)
                disp("skipped region without ephys: "+string(idmap.ccfid2reg(sink_ccfid(ii))));
                continue
            end
            allen_mat=[allen_mat;allen_density];
            allen_meta=[allen_meta;1,ii,sink_ccfid(ii)];
        end
    end
    % afferent_proj_dense(sink) ~ feature_proportion
    idx4corr=cell2mat(sink_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
    for ii=1:numel(src_ccfid)
        allen_density=sink_src_mat(idx4corr,ii); % from one alternating target, to idx4corr
        if nnz(allen_density~=0)<8
            disp(ii)
            continue
        elseif ~ismember(idmap.ccfid2reg(src_ccfid(ii)),map_regs)
            disp("skipped region without ephys: "+string(idmap.ccfid2reg(src_ccfid(ii))));
            continue
        end
        
        allen_mat=[allen_mat;allen_density.'];
        allen_meta=[allen_meta;2,ii,src_ccfid(ii)];
    end

    intersect_ccfid=cell2mat(idmap.reg2ccfid.values(intersect_regs));

    %% one region corr, bar plot
    
    if opt.corr1
        
        for regii=1:size(allen_mat,1)
            % TODO: remove injection site
            switch(corr_type)
                case 'Spearman'
                    vsel=intersect_ccfid~=allen_meta(regii,3);
                    [r,p]=corr(reshape(allen_mat(regii,vsel).',[],1),...
                        reshape(feat_prop(vsel),[],1),'type','Spearman');
                case 'PearsonLogLog'
                    vsel=(allen_mat(regii,:).')>0 & feat_prop>0 & intersect_ccfid~=allen_meta(regii,3);
                    [r,p]=corr(reshape(log10(allen_mat(regii,vsel)),[],1),...
                        reshape(log10(feat_prop(vsel)),[],1),'type','Pearson');
                case 'PearsonLinearLog'
                    vsel=(allen_mat(regii,:).')>0 & intersect_ccfid~=allen_meta(regii,3);
                    [r,p]=corr(reshape(log10(allen_mat(regii,vsel)),[],1),...
                        reshape(feat_prop(vsel),[],1),'type','Pearson');
            end
            one_reg_corr_list=[one_reg_corr_list;featii,epochii,regii,double(allen_meta(regii,:)),r,p];
            %====================================^^^1^^^^^^2^^^^^^3^^^^^^^^^^^4,5,6^^^^^^^^^^^^^7^8^^
        end
%         save('one_reg_corr_list.mat','one_reg_corr_list')

        if opt.plot1 % plot one-region: coding proportion correlation bars
            s_list=sortrows(one_reg_corr_list(one_reg_corr_list(:,1)==featii & one_reg_corr_list(:,2)==epochii,:),7,'ascend');
%             s_list=s_list([1:2,end-1:end],:);
           
            xlbl=idmap.ccfid2reg.values(num2cell(s_list(:,6)));
            xlbl=cellfun(@(x) x{1}, xlbl,'UniformOutput',false);

            fh.(sprintf('feat%d',featii))=figure('Color','w','Position',[32,32,800,800]);
            tiledlayout(4,2);

            nexttile(5,[1,2])
            hold on
            %                 bhto=bar(find(s_list(:,4)==1),s_list(s_list(:,4)==1,7),0.6,'white');
            %                 bhfrom=bar(find(s_list(:,4)==2),s_list(s_list(:,4)==2,7),0.6,'black');
            bh=bar(s_list(:,7),0.8,'FaceColor','flat','Horizontal','off');
            if false
                if false
                    fid=fopen(fullfile('binary','upload','F1N_fraction_anatomy_correlation.json'),'w');
                elseif 0>1
                    fid=fopen(fullfile('binary','upload','F1SBelow_region_timing_anatomy_correlation.json'),'w');
                else
                    fid=fopen(fullfile('binary','upload','SF4B_region_timing_anatomy_correlation.json'),'w');
                end
                fprintf(fid,jsonencode(table(xlbl,s_list(:,7),'VariableNames',{'Region','Pearsons_r'})));
                fclose(fid);
            end


            bh.CData(s_list(:,4)==1,:)=repmat([0.5,0.5,0.5],nnz(s_list(:,4)==1),1);
            bh.CData(s_list(:,4)==2,:)=repmat([0,0,0],nnz(s_list(:,4)==2),1);
            ylabel('Connectivity-coding proportion correlation (Pearson''s r)')
            set(gca(),'XTick',1:size(s_list,1),'XTickLabel',xlbl,'XTickLabelRotation',90);
            
            sgtitle(fns{featii},'Interpreter','none');
%                 title(sprintf('epoch%d-feature%d',epochii,featii))
            ylim([-1,1]);

            for rr=1:numel(xlbl)
                c=ephys.getRegColor(xlbl{rr},'large_area',true);
                text(rr,1,xlbl{rr},'HorizontalAlignment','center','VerticalAlignment','bottom','Color',c,'Rotation',90)
            end            
            
            %   legend([bhto,bhfrom],{'Connectivity to this region','Connectivity from this region'})
            %   exportgraphics(fh.bar1,sprintf('frac_allen_mdls_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')

%             nexttile(1,[1,3])
            nexttile(1,[2,1]) % positive max
            hold on
            glmidx=find(allen_meta(:,1)==s_list(1,4) & allen_meta(:,2)==s_list(1,5));
            rmv_inj=intersect_ccfid~=allen_meta(s_list(1,3),3);
            for ll=1:numel(feat_prop)
                if intersect_ccfid(ll)==allen_meta(s_list(1,3),3)
                    continue
                end
                c=ephys.getRegColor(intersect_regs{ll},'large_area',true);
                scatter(allen_mat(glmidx,ll),feat_prop(ll),4,c,'filled','o')
                text(allen_mat(glmidx,ll),feat_prop(ll),intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
            end

            
            xlabel('Projection density (px/px)');
            ylabel('Regional averaged feature index')
            if strcmp(corr_type,'PearsonLinearLog')
                set(gca,'XScale','log','YScale','linear')
                coord=[log10(allen_mat(glmidx,rmv_inj).'),feat_prop(rmv_inj)];
                coord(:,3)=1;
                regres=coord(:,[1,3])\coord(:,2);
                xx=minmax(allen_mat(glmidx,rmv_inj));
                yy=log10(xx).*regres(1)+regres(2);
                plot(xx,yy,'--k');
            else
                set(gca,'XScale','log','YScale','log')
                coord=log10([allen_mat(glmidx,rmv_inj).',feat_prop(rmv_inj)]);
                coord(:,3)=1;
                regres=coord(:,[1,3])\coord(:,2);
                xx=minmax(allen_mat(glmidx,rmv_inj));
                yy=10.^(log10(xx).*regres(1)+regres(2));
                plot(xx,yy,'--k');
            end
            title(sprintf(' r = %.3f, p = %.3f',s_list(1,7),s_list(1,8)));
            if false
                if false
                    fid=fopen(fullfile('binary','upload','F1M_fraction_AON_projection_correlation.json'),'w');
                    ttbl=table(intersect_regs,allen_mat(glmidx,:).',feat_prop,'VariableNames',{'Region','Projection_density','WM_neuron_proportion'});
                else
                    fid=fopen(fullfile('binary','upload','F1SAbove_region_wave_timing_RHP_projection_correlation.json'),'w');
                    ttbl=table(intersect_regs,allen_mat(glmidx,:).',feat_prop,'VariableNames',{'Region','Projection_density','region_wave_timing'});
                end
                fprintf(fid,jsonencode(ttbl));
                fclose(fid);
            end            

            
            nexttile(2,[2,1]) %negative max
            hold on
            glmidx=find(allen_meta(:,1)==s_list(end,4) & allen_meta(:,2)==s_list(end,5));
            rmv_inj=intersect_ccfid~=allen_meta(s_list(end,3),3);
            for ll=1:numel(feat_prop)
                if intersect_ccfid(ll)==allen_meta(s_list(end,3),3)
                    continue
                end
                c=ephys.getRegColor(intersect_regs{ll},'large_area',true);
                scatter(allen_mat(glmidx,ll),feat_prop(ll),4,c,'filled','o')
                text(allen_mat(glmidx,ll),feat_prop(ll),intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
            end
            xlabel('Projection density (px/px)');
            ylabel('Regional averaged feature index')
            if strcmp(corr_type,'PearsonLinearLog')
                set(gca,'XScale','log','YScale','linear')
                coord=[log10(allen_mat(glmidx,rmv_inj).'),feat_prop(rmv_inj)];
                coord(:,3)=1;
                regres=coord(:,[1,3])\coord(:,2);
                xx=minmax(allen_mat(glmidx,rmv_inj));
                yy=log10(xx).*regres(1)+regres(2);
                plot(xx,yy,'--k');
            else
                set(gca,'XScale','log','YScale','log')
                coord=log10([allen_mat(glmidx,rmv_inj).',feat_prop(rmv_inj)]);
                coord(:,3)=1;
                regres=coord(:,[1,3])\coord(:,2);
                xx=minmax(allen_mat(glmidx,rmv_inj));
                yy=10.^(log10(xx).*regres(1)+regres(2));
                plot(xx,yy,'--k');
            end
            title(sprintf(' r = %.3f, p = %.3f',s_list(end,7),s_list(end,8)));
            if false
                fid=fopen(fullfile('binary','upload','F1M_fraction_AON_projection_correlation.json'),'w');
                ttbl=table(intersect_regs,allen_mat(glmidx,:).',feat_prop,'VariableNames',{'Region','Projection_density','WM_neuron_proportion'});
                fprintf(fid,jsonencode(ttbl));
                fclose(fid);
            end

            th=nexttile(7,[1,2]); % data table
            tbl=cell(0);
            if isfield(opt,'stats_type') && ~isempty(opt.stats_type)
                tbl=[tbl;'Stats';opt.stats_type];
            end
            if isfield(opt,'data_type') && ~isempty(opt.data_type)
                tbl=[tbl;'Data';opt.data_type];
            end
            tbl=[tbl;'Range';opt.range];
            ephys.util.figtable(gcf(),th,tbl)

%             exportgraphics(fh.corr1,sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')

        end
    end
    %% two region glm interaction
    if opt.corr2
        comb2=nchoosek(1:size(allen_mat,1),2);
        mdlid_rsq_AICC=[];
        for jj=1:size(comb2,1)
            if rem(jj,1000)==0
                disp(jj)
            end
            allen_A=log10(allen_mat(comb2(jj,1),:).'); % from one alternating target, to idx4corr
            allen_B=log10(allen_mat(comb2(jj,2),:).');
            mdl=fitglm([allen_A,allen_B],feat_prop,'linear'); % (1|2)jj => glmxmat==glmxmeta => (sink_ccfid|src_ccfid)
%             two_reg_corr_list=[two_reg_corr_list;featii,epochii,jj,double(glmxmeta(comb2(jj,1),:)),double(glmxmeta(comb2(jj,2),:)),mdl.Coefficients.Estimate.',mdl.Rsquared.Ordinary,mdl.ModelCriterion.AICc];
            %============================================^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            mdlid_rsq_AICC=[mdlid_rsq_AICC;jj,mdl.Rsquared.Ordinary,mdl.ModelCriterion.AICc,mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(3)];
        end
        if opt.plot2
            mdlid_rsq_AICC=sortrows(mdlid_rsq_AICC,3);
            signs=["-","+","+"];
            ytk=cell(0,0);
            for kk=1:10
                maxidx=mdlid_rsq_AICC(kk,1);
%                 if glmxmeta(comb2(maxidx,1),1)==1
%                     reg1=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,1),2)));
%                     reg1=['To ',reg1{1}];
%                 else
                    reg1=(idmap.ccfid2reg(src_ccfid(allen_meta(comb2(maxidx,1),2))));
                    reg1=signs(sign(mdlid_rsq_AICC(kk,4))+2)+string(reg1);
%                 end
%                 if glmxmeta(comb2(maxidx,2),1)==1
%                     reg2=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,2),2)));
%                     reg2=['To ',reg2{1}];
%                 else
                    reg2=idmap.ccfid2reg(src_ccfid(allen_meta(comb2(maxidx,2),2)));
                    reg2=signs(sign(mdlid_rsq_AICC(kk,5))+2)+string(reg2);
%                 end
%                 disp([reg1,',',reg2,' ',num2str(mdlid_rsq_AICC(kk,2:3))])
                ytk{end+1}=reg1+", "+reg2;
            end

            fh.(sprintf('feat%d_2f',featii))=figure('Color','w','Position',[32,32,400,600]);
            tiledlayout(3,1)
            nexttile(3);
            bar(sqrt(mdlid_rsq_AICC(1:10,2)),'Horizontal','on','FaceColor','k')
            set(gca(),'YDir','reverse','YTick',1:10,'YTickLabel',ytk(1:10))
            xlabel('Person''s r')
            xlim([0,1])
            title(sprintf('selectivity %d - epoch %d',featii,epochii));
%             exportgraphics(fh,sprintf('frac_allen_mdls_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')
%             close(fh)

            maxidx=mdlid_rsq_AICC(1,1);
            allen_A=log10(allen_mat(comb2(maxidx,1),:).'); % from one alternating target, to idx4corr
            allen_B=log10(allen_mat(comb2(maxidx,2),:).');
            mdl=fitglm([allen_A,allen_B],feat_prop,'linear');


%             disp(sqrt(mdl.Rsquared.Ordinary));
%             if glmxmeta(comb2(maxidx,1),1)==1
%                 reg1=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,1),2)));
%                 reg1=['To_',reg1{1}];
%             else
                reg1=(idmap.ccfid2reg(src_ccfid(allen_meta(comb2(maxidx,1),2))));
                reg1=['From_',reg1{1}];
%             end
%             if glmxmeta(comb2(maxidx,2),1)==1
%                 reg2=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(maxidx,2),2)));
%                 reg2=['To_',reg2{1}];
%             else
                reg2=idmap.ccfid2reg(src_ccfid(allen_meta(comb2(maxidx,2),2)));
                reg2=['From_',reg2{1}];
%             end
            nexttile(1,[2,1])
%             fh=figure('Color','w','Position',[32,32,320,320]);
            hold on
            for ll=1:numel(feat_prop)
                c=ephys.getRegColor(intersect_regs{ll},'large_area',true);
                scatter(mdl.Fitted.Response(ll),feat_prop(ll),4,c,'filled','o')
                text(mdl.Fitted.Response(ll),feat_prop(ll),intersect_regs{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
            end
%             title(sprintf('Coding proportion ~ %.3f [%s] + %.3f [%s] +%.3f [%s]', ...
%                 mdl.Coefficients.Estimate(2),reg1,mdl.Coefficients.Estimate(3),reg2,mdl.Coefficients.Estimate(4),'interaction'), ...
            title(sprintf('Coding proportion ~ %.3f [%s] + %.3f [%s]', ...
                mdl.Coefficients.Estimate(2),reg1,mdl.Coefficients.Estimate(3),reg2), ...
            'Interpreter','none')
            xlabel(sprintf('Model %d prediction',maxidx))
            ylabel('Dependent feature')
            if strcmp(corr_type,'PearsonLinearLog')
                set(gca,'XScale','log','YScale','linear')
            else
                set(gca,'XScale','log','YScale','log')
            end

            text(min(xlim()),max(ylim()),sprintf(' r = %.3f, AIC = %.1f',sqrt(mdl.Rsquared.Ordinary),mdl.ModelCriterion.AIC),'HorizontalAlignment','left','VerticalAlignment','top');
            % xlim([0.15,0.5])
            % ylim([0.15,0.5])
%             exportgraphics(fh,sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')
%             close(fh)
            %         keyboard()

        end
%         save('two_reg_corr_list.mat','two_reg_corr_list')
    end
    %%
    if opt.corr3
        comb3=nchoosek(1:size(allen_mat,1),3);
        max3rsqp=0;
        max3idxp=-1;
        for jj=127706%1:size(comb3,1)
            if rem(jj,1000)==0
                disp(jj)
            end
            allen_A=allen_mat(comb3(jj,1),:).'; % from one alternating target, to idx4corr
            allen_B=allen_mat(comb3(jj,2),:).';
            allen_C=allen_mat(comb3(jj,3),:).';
            mdl=fitglm([allen_A,allen_B,allen_C],feat_prop,'poly111');
            if mdl.Rsquared.Ordinary>max3rsqp
                max3rsqp=mdl.Rsquared.Ordinary;
                max3idxp=jj;
            end

            if mdl.Rsquared.Ordinary>0.64
                disp(jj)
                disp(sqrt(mdl.Rsquared.Ordinary));
                if allen_meta(comb3(jj,1),1)==1
                    reg1=idmap.ccfid2reg(sink_ccfid(allen_meta(comb3(jj,1),2)));
                    reg1=['To_',reg1{1}];
                else
                    reg1=(idmap.ccfid2reg(src_ccfid(allen_meta(comb3(jj,1),2))));
                    reg1=['From_',reg1{1}];
                end
                if allen_meta(comb3(jj,2),1)==1
                    reg2=idmap.ccfid2reg(sink_ccfid(allen_meta(comb3(jj,2),2)));
                    reg2=['To_',reg2{1}];
                else
                    reg2=idmap.ccfid2reg(src_ccfid(allen_meta(comb3(jj,2),2)));
                    reg2=['From_',reg2{1}];
                end

                if allen_meta(comb3(jj,3),1)==1
                    reg3=idmap.ccfid2reg(sink_ccfid(allen_meta(comb3(jj,3),2)));
                    reg3=['To_',reg3{1}];
                else
                    reg3=idmap.ccfid2reg(src_ccfid(allen_meta(comb3(jj,3),2)));
                    reg3=['From_',reg3{1}];
                end
                if opt.plot3
                    fh=figure('Color','w','Position',[32,32,900,650]);
                    scatter(mdl.Fitted.Response,feat_prop,4,'red','filled','o')
                    text(mdl.Fitted.Response,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
                    title(sprintf('Coding proportion ~ %.3f [%s] + %.3f [%s] +%.3f [%s]', ...
                        mdl.Coefficients.Estimate(2),reg1,mdl.Coefficients.Estimate(3),reg2,mdl.Coefficients.Estimate(4),reg3), ...
                        'Interpreter','none')
                    xlabel(sprintf('Model %d prediction',jj))
                    ylabel('Proportion of coding neuron')
                    text(min(xlim()),max(ylim()),sprintf(' r = %.3f, AIC = %.1f',sqrt(mdl.Rsquared.Ordinary),mdl.ModelCriterion.AIC),'HorizontalAlignment','left','VerticalAlignment','top');
                end    %         keyboard()
            end
        end
    end
end