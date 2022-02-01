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
    trlCount=[];
    for skey=reshape(fieldnames(ordmat),1,[])
        trls=ordmat.(skey{1}).block_meta;
        csel=trls(:,9)~=0 & trls(:,10)~=0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,11),1:4);
        sess_lbl=arrayfun(@(x) sprintf('b%dt%d',trls(x,8),trls(x,11)),find(csel),'UniformOutput',false);
        frs=mean(ordmat.(skey{1}).fr(csel,:,1:10),3); % [nTrl,nSU]
        trlTypes=["b3t1","b3t2","b3t3","b3t4","b6t1","b6t2","b6t3","b6t4"];
        trlSess=nan(1,8);
        for trlIdx=1:8
            trlType=trlTypes(trlIdx);
            trlSess(trlIdx)=nnz(strcmp(sess_lbl,trlType));
            if ~isfield(frmap,trlType)
                frmap.(trlType)=containers.Map('KeyType','char','ValueType','any');
            end
            for suidx=1:size(frs,2) % TODO: proper batch init
                frmap.(trlType)(sprintf('%su%d',skey{1},suidx))=frs(strcmp(sess_lbl,trlType),suidx);
            end
        end
        trlCount=[trlCount;trlSess];
    end

end
%% available trials for least repeated condition
%figure();histogram(trlCount(:),1:32)


%% circle of PCA
if true
    sukeys=frmap.b3t1.keys();
    trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),frmap.b3t1.keys()),find(min(trlCount,[],2)>rpt));
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
rpt=10;
sukeys=frmap.b3t1.keys();
trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),frmap.b3t1.keys()),find(min(trlCount,[],2)>rpt));
sukeys=sukeys(trlSel);

cv=cvpartition(rpt,'KFold',10);
dmat=struct();
for lbl=["b3t1","b3t2","b3t3","b3t4","b6t1","b6t2","b6t3","b6t4"]
    dmat.(lbl)=cell2mat(cellfun(@(x) randsample(frmap.(lbl)(x),rpt),sukeys,'UniformOutput',false));
end
[result,shuf]=deal([]);
for kf=1:cv.NumTestSets
    [xx_train,yy_train,xx_test,yy_test]=deal([]);
    for lbl=["b3t1","b3t2","b3t3","b3t4","b6t1","b6t2","b6t3","b6t4"]
        %         keyboard();
        xx_train=[xx_train;dmat.(lbl)(training(cv,kf),:)];
        xx_test=[xx_test;dmat.(lbl)(test(cv,kf),:)];
        yy_train=[yy_train;repmat(lbl,nnz(training(cv,kf)),1)];
        yy_test=[yy_test;repmat(lbl,nnz(test(cv,kf)),1)];

    end
    norm_train=normalize(xx_train,'center',C,'scale',S);
    pc_train=norm_train*coeff;
    SVMM=fitcecoc(pc_train,yy_train);
    norm_test=normalize(xx_test,'center',C,'scale',S);
    pc_test=norm_test*coeff;
    modelPredict=SVMM.predict(pc_test);
    result=[result;modelPredict==yy_test];
    shuf=[shuf;randsample(modelPredict,numel(modelPredict))==yy_test];
end
