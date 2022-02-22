%% efferent_proj_dense(src) ~ feature_proportion

allen_src_regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(src_ccfid)),'UniformOutput',false);
intersect_regs=intersect(allen_src_regs,chregs);

idx4corr=cell2mat(src_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
feat_prop_cell=feat_reg_map.values(intersect_regs);
feat_prop=cellfun(@(x) x(1),feat_prop_cell);

[glmxmat,glmxmeta]=deal([]);
for ii=1:numel(sink_ccfid)
    allen_density=log10(sink_src_mat(ii,idx4corr)+1e-12); % from idx4corr, to one alternating target
    glmxmat=[glmxmat;allen_density];
    glmxmeta=[glmxmeta;1,ii];
end

%% afferent_proj_dense(sink) ~ feature_proportion
idx4corr=cell2mat(sink_idx_map.values(idmap.reg2ccfid.values(intersect_regs)));
for ii=1:numel(src_ccfid)
    allen_density=log10(sink_src_mat(idx4corr,ii)+1e-12); % from one alternating target, to idx4corr
    glmxmat=[glmxmat;allen_density.'];
    glmxmeta=[glmxmeta;2,ii];
end
%%

comb2=nchoosek(1:size(glmxmat,1),2);
maxrsq=0;
maxidx=-1;
for jj=3992%:size(comb2,1)
    if rem(jj,1000)==0
        disp(jj)
    end
    allen_A=glmxmat(comb2(jj,1),:).'; % from one alternating target, to idx4corr
    allen_B=glmxmat(comb2(jj,2),:).';
    mdl=fitglm([allen_A,allen_B],feat_prop,'interactions');
    if mdl.Rsquared.Ordinary>maxrsq
        maxrsq=mdl.Rsquared.Ordinary;
        maxidx=jj;
    end
%     if mdl.Rsquared.Ordinary>0.64
        disp(sqrt(mdl.Rsquared.Ordinary));
        if glmxmeta(comb2(jj,1),1)==1
            reg1=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(jj,1),2)));
            reg1=['To_',reg1{1}];
        else
            reg1=(idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(jj,1),2))));
            reg1=['From_',reg1{1}];
        end
        if glmxmeta(comb2(jj,2),1)==1
            reg2=idmap.ccfid2reg(sink_ccfid(glmxmeta(comb2(jj,2),2)));
            reg2=['To_',reg2{1}];
        else
            reg2=idmap.ccfid2reg(src_ccfid(glmxmeta(comb2(jj,2),2)));
            reg2=['From_',reg2{1}];
        end

        fh=figure('Color','w','Position',[32,32,900,650]);
        scatter(mdl.Fitted.Response,feat_prop,4,'red','filled','o')
        text(mdl.Fitted.Response,feat_prop,intersect_regs,'HorizontalAlignment','center','VerticalAlignment','top');
        title(sprintf('Coding proportion ~ %.3f [%s] + %.3f [%s] +%.3f [%s]', ...
            mdl.Coefficients.Estimate(2),reg1,mdl.Coefficients.Estimate(3),reg2,mdl.Coefficients.Estimate(4),'interaction'), ...
            'Interpreter','none')
        xlabel(sprintf('Model %d prediction',jj))
        ylabel('Proportion of coding neuron')
        text(min(xlim()),max(ylim()),sprintf(' r = %.3f, AIC = %.1f',sqrt(mdl.Rsquared.Ordinary),mdl.ModelCriterion.AIC),'HorizontalAlignment','left','VerticalAlignment','top');
        keyboard()
%     end
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

%     if mdl.Rsquared.Ordinary>0.64
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
%     end
end

