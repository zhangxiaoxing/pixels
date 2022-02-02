%TODO brain region filter, olfaction filter.
% function out=get_ordinal_mat(opt)
%% gen data
if false
    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta();
    homedir=ephys.util.getHomedir('type','raw');
    for ii=reshape(cell2mat(sessmap.keys()),1,[])
        disp(ii)
        fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        dur_resp=behav.tag_block(trials,'wt',false);
        block_meta=[trials,dur_resp(:,end)];
        sesssel=meta.sess==ii;
        regsel=strcmp(meta.reg_tree(2,sesssel),'CTX');

        fr=fr(:,regsel,:);

        out.(sprintf('s%d',ii)).fr=fr;
        out.(sprintf('s%d',ii)).block_meta=block_meta;
    end
    ordmat=out;
    frmap=struct();
    [c_trlCount,e_trlCount]=deal([]);
    for skey=reshape(fieldnames(ordmat),1,[])
        trls=ordmat.(skey{1}).block_meta;
        csel=trls(:,9)~=0 & trls(:,10)~=0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,11),1:4);
        esel=trls(:,10)==0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,11),1:4);
        c_sess_lbl=arrayfun(@(x) sprintf('b%dt%d',trls(x,8),trls(x,11)),find(csel),'UniformOutput',false);
        e_sess_lbl=arrayfun(@(x) sprintf('b%dt%d',trls(x,8),trls(x,11)),find(esel),'UniformOutput',false);
        c_frs=mean(ordmat.(skey{1}).fr(csel,:,1:10),3); % [nTrl,nSU]
        e_frs=mean(ordmat.(skey{1}).fr(esel,:,1:10),3); % [nTrl,nSU]
        trlTypes=["b3t1","b3t2","b3t3","b3t4","b6t1","b6t2","b6t3","b6t4"];
        [c_trlSess,e_trlSess]=deal(nan(1,8));
        for trlIdx=1:8
            trlType=trlTypes(trlIdx);
            c_trlSess(trlIdx)=nnz(strcmp(c_sess_lbl,trlType));
            e_trlSess(trlIdx)=nnz(strcmp(e_sess_lbl,trlType));
            if ~isfield(frmap,trlType)
                frmap.(trlType)=containers.Map('KeyType','char','ValueType','any');
                frmap.(strjoin([trlType,"_e"],''))=containers.Map('KeyType','char','ValueType','any');
            end
            for suidx=1:size(c_frs,2) % TODO: vectorized batch init
                frmap.(trlType)(sprintf('%su%d',skey{1},suidx))=c_frs(strcmp(c_sess_lbl,trlType),suidx);
                frmap.(strjoin([trlType,"_e"],''))(sprintf('%su%d',skey{1},suidx))=e_frs(strcmp(e_sess_lbl,trlType),suidx);
            end
        end
        c_trlCount=[c_trlCount;c_trlSess];
        e_trlCount=[e_trlCount;e_trlSess];
    end

end
%% available trials for least repeated condition
%figure();histogram(trlCount(:),1:32)


%% circle of PCA
min_trl=10;
min_e_trl=1;
if false
    sukeys=frmap.b3t1.keys();
    trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),frmap.b3t1.keys()),find(min(c_trlCount,[],2)>min_trl & min(e_trlCount,[],2)>=min_e_trl));
    sukeys=sukeys(trlSel);
    pcamat=[...
        cellfun(@(x) mean(x),frmap.b3t1.values(sukeys));...
        cellfun(@(x) mean(x),frmap.b3t2.values(sukeys));...
        cellfun(@(x) mean(x),frmap.b3t3.values(sukeys));...
        cellfun(@(x) mean(x),frmap.b3t4.values(sukeys));...
        cellfun(@(x) mean(x),frmap.b6t1.values(sukeys));...
        cellfun(@(x) mean(x),frmap.b6t2.values(sukeys));...
        cellfun(@(x) mean(x),frmap.b6t3.values(sukeys));...
        cellfun(@(x) mean(x),frmap.b6t4.values(sukeys))];

    [normmat,C,S]=normalize(pcamat,1,"range");

    [coeff,score,latent]=pca(normmat);
    figure()
    for ii=1:7
        subplot(2,4,ii)
        plot(score(:,ii))
    end

    fh=figure('Color','w','Position',[100,100,250,225]);
    hold on
    plot(score([1:8,1],1),score([1:8,1],2),'-','Color',[0.5,0.5,0.5])

    ph3=plot(score(1:4,1),score(1:4,2),'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerSize',9);
    ph6=plot(score(5:8,1),score(5:8,2),'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',9);

    tidx=[1:4,1:4];
    arrayfun(@(x) text(score(x,1),score(x,2),sprintf('No.%d',tidx(x)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10),1:8)
    xlabel('PC1')
    ylabel('PC2')
    legend([ph3,ph6],{'3s delay','6s delay'},'Location','northoutside','Orientation','horizontal')
    xlim([-30,25])
    exportgraphics(fh,'Block_trial_PC.pdf','ContentType','vector')
end
%% SVM

sukeys=frmap.b3t1.keys();
trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),frmap.b3t1.keys()),find(min(c_trlCount,[],2)>min_trl & min(e_trlCount,[],2)>=min_e_trl));
sukeys=sukeys(trlSel);

cv=cvpartition(min_trl,'KFold',10);
dmat=struct();
for trlType=["b3t1","b3t2","b3t3","b3t4","b6t1","b6t2","b6t3","b6t4"]
    dmat.(trlType)=cell2mat(cellfun(@(x) datasample(frmap.(trlType)(x),min_trl),sukeys,'UniformOutput',false));
    dmat.(strjoin([trlType,"_e"],''))=cell2mat(cellfun(@(x) datasample(frmap.(strjoin([trlType,"_e"],''))(x),min_e_trl),sukeys,'UniformOutput',false));
end
[result,shuf,result_e]=deal([]);
for kf=1:cv.NumTestSets
    [xx_train,yy_train,xx_test,yy_test,xx_e_test,yy_e_test]=deal([]);
    for trlType=["b3t1","b3t2","b3t3","b3t4","b6t1","b6t2","b6t3","b6t4"]
        %         keyboard();
        xx_train=[xx_train;dmat.(trlType)(training(cv,kf),:)];
        xx_test=[xx_test;dmat.(trlType)(test(cv,kf),:)];
        xx_e_test=[xx_e_test;dmat.(strjoin([trlType,"_e"],''))];
        yy_train=[yy_train;repmat(trlType,nnz(training(cv,kf)),1)];
        yy_test=[yy_test;repmat(trlType,nnz(test(cv,kf)),1)];
        yy_e_test=[yy_e_test;repmat(trlType,min_e_trl,1)];
    end
    norm_train=normalize(xx_train,'center',C,'scale',S);
    pc_train=norm_train*coeff;
    SVMM=fitcecoc(pc_train,yy_train);
    norm_test=normalize(xx_test,'center',C,'scale',S);
    pc_test=norm_test*coeff;
    modelPredict=SVMM.predict(pc_test);
    result=[result;modelPredict==yy_test];
    shuf=[shuf;randsample(modelPredict,numel(modelPredict))==yy_test];

    norm_e_test=normalize(xx_e_test,'center',C,'scale',S);
    pc_e=norm_e_test*coeff;
    e_predict=SVMM.predict(pc_e);
    result_e=[result_e;e_predict==yy_e_test];
end
