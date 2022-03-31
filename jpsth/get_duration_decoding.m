%TODO brain region filter, olfaction filter.
function out=get_duration_decoding(opt)
arguments
    opt.new_data (1,1) logical=true
%     opt.plot_PCA (1,1) logical=false
    opt.calc_dec (1,1) logical=true
    opt.plot_dec (1,1) logical=true
end
%% gen data

if opt.new_data
    dur_mat=struct();
    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta();
    homedir=ephys.util.getHomedir('type','raw');
    for ii=reshape(cell2mat(sessmap.keys()),1,[])
        disp(ii)
        fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
%         suid=h5read(fpath,'/SU_id');

        sesssel=meta.sess==ii;
        regsel=ismember(meta.reg_tree(1,sesssel),{'CH','BS'});

        fr=fr(:,regsel,:);

        fr_t_align=nan(size(fr,1),size(fr,2)); %trial, su, pre-test bin-averaged
        fr_t_align(trials(:,8)==3,:)=mean(fr(trials(:,8)==3,:,17:28),3); % early/late delay
        fr_t_align(trials(:,8)==6,:)=mean(fr(trials(:,8)==6,:,17:28),3); % early/late delay
        

        dur_mat.(sprintf('s%d',ii)).fr_t=fr_t_align;
        dur_mat.(sprintf('s%d',ii)).trials=trials;
    end

    frmap=struct();
    frmap.d3=containers.Map('KeyType','char','ValueType','any');
    frmap.d6=containers.Map('KeyType','char','ValueType','any');
    frmap.d3_e=containers.Map('KeyType','char','ValueType','any');
    frmap.d6_e=containers.Map('KeyType','char','ValueType','any');

    [c_trlCount,e_trlCount]=deal([]);
    for skey=reshape(fieldnames(dur_mat),1,[])
        trls=dur_mat.(skey{1}).trials;
        csel=trls(:,9)~=0 & trls(:,10)~=0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        esel=trls(:,10)==0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        c_dur_lbl=trls(csel,8);
        e_dur_lbl=trls(esel,8);
        c_frs=dur_mat.(skey{1}).fr_t(csel,:); % [nTrl,nSU]
        e_frs=dur_mat.(skey{1}).fr_t(esel,:); % [nTrl,nSU]
        [c_trlSess,e_trlSess]=deal([0,0]);
        for trl_dur=[3,6]
            c_trlSess(trl_dur/3)=nnz(c_dur_lbl==trl_dur);
            e_trlSess(trl_dur/3)=nnz(e_dur_lbl==trl_dur);
            for suidx=1:size(c_frs,2) % TODO: vectorized batch init
                frmap.(sprintf('d%d',trl_dur))(sprintf('%su%d',skey{1},suidx))=c_frs(c_dur_lbl==trl_dur,suidx);
                frmap.(sprintf('d%d_e',trl_dur))(sprintf('%su%d',skey{1},suidx))=e_frs(e_dur_lbl==trl_dur,suidx);
            end
        end
        c_trlCount=[c_trlCount;c_trlSess];
        e_trlCount=[e_trlCount;e_trlSess];
    end
    save('dur_fr.mat','frmap','c_trlCount','e_trlCount');
else
    load('dur_fr.mat','frmap','c_trlCount','e_trlCount');
end
%% available trials for least repeated condition
%figure();histogram(trlCount(:),1:32)


%% decoding

min_trl=10;
min_e_trl=1;
sukeys_all=frmap.d3.keys();
trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),frmap.d3.keys()),find(min(c_trlCount,[],2)>min_trl & min(e_trlCount,[],2)>=min_e_trl));
sukeys_all=sukeys_all(trlSel);


if opt.calc_dec
    out=struct();
    for n_su=[50,100,500,1000]
        [result,shuf,result_e]=deal([]);
        for resamp_rpt=1:50%15
            sukeys=datasample(sukeys_all,n_su);
            pcamat=[...
                cellfun(@(x) min(x),frmap.d3.values(sukeys));...
                cellfun(@(x) max(x),frmap.d3.values(sukeys));...
                cellfun(@(x) min(x),frmap.d6.values(sukeys));...
                cellfun(@(x) max(x),frmap.d6.values(sukeys))];

             [normmat,C,S]=normalize(pcamat,1,"range");
%              [coeff,~,~]=pca(normmat);

            cv=cvpartition(min_trl,'KFold',10);
            dmat=struct();
            for trlType=["d3","d6"]
                dmat.(trlType)=cell2mat(cellfun(@(x) datasample(frmap.(trlType)(x),min_trl),sukeys,'UniformOutput',false));
                dmat.(strjoin([trlType,"_e"],''))=cell2mat(cellfun(@(x) datasample(frmap.(strjoin([trlType,"_e"],''))(x),min_e_trl),sukeys,'UniformOutput',false));
            end

            for kf=1:cv.NumTestSets
                [xx_train,yy_train,xx_test,yy_test,xx_e_test,yy_e_test]=deal([]);
                for trlType=["d3","d6"]
                    xx_train=[xx_train;dmat.(trlType)(training(cv,kf),:)];
                    xx_test=[xx_test;dmat.(trlType)(test(cv,kf),:)];
                    xx_e_test=[xx_e_test;dmat.(strjoin([trlType,"_e"],''))];
                    yy_train=[yy_train;repmat(trlType,nnz(training(cv,kf)),1)];
                    yy_test=[yy_test;repmat(trlType,nnz(test(cv,kf)),1)];
                    yy_e_test=[yy_e_test;repmat(trlType,min_e_trl,1)];
                end
                norm_train=normalize(xx_train,'center',C,'scale',S);
%                 pc_train=norm_train*coeff;
                SVMM=fitcecoc(norm_train,yy_train);
                norm_test=normalize(xx_test,'center',C,'scale',S);
%                 pc_test=norm_test*coeff;
                modelPredict=SVMM.predict(norm_test);
                result=[result;modelPredict==yy_test];
                shuf=[shuf;randsample(modelPredict,numel(modelPredict))==yy_test];

                norm_e_test=normalize(xx_e_test,'center',C,'scale',S);
%                 pc_e=norm_e_test*coeff;
                e_predict=SVMM.predict(norm_e_test);
                result_e=[result_e;e_predict==yy_e_test];
            end
        end
        out.(sprintf('result_%dsu',n_su))=result;
        out.(sprintf('shuf_%dsu',n_su))=shuf;
        out.(sprintf('etrial_%dsu',n_su))=result_e;
    end
    save('duration_decoding.mat','out');
elseif opt.plot_dec
    load('duration_decoding.mat','out');
end

if opt.plot_dec
    colors=containers.Map(["result","shuf","etrial"],{'r','k','b'});
    fh=figure('Color','w','Position',[100,100,250,225]);
    hold on;
    for ptype=["result","shuf","etrial"]
        %     datamat=cell2mat(arrayfun(@(x) out.(sprintf('%s_%dsu',ptype,x)),[50,100,200,500,1000],'UniformOutput',false));
        n_su=[50,100,500,1000];
        phat=nan(1,numel(n_su));
        pci=nan(2,numel(n_su));
        for nidx=1:numel(n_su)
            dd=out.(sprintf('%s_%dsu',ptype,n_su(nidx)));
            [phat(nidx),pci(:,nidx)]=binofit(nnz(dd),numel(dd));
        end
        fill([n_su,fliplr(n_su)],[pci(1,:),fliplr(pci(2,:))],colors(ptype),'EdgeColor','none','FaceAlpha',0.1);
        plot(n_su,phat,'Color',colors(ptype))
    end
    xlabel('Number of neurons')
    ylabel('Classification accuracy (%)')
    ylim([0,1])
    set(gca(),'YTick',0:0.1:1,'YTickLabel',0:10:100)
    xlim([0,max(n_su)])
    exportgraphics(fh,'duration_decoding.pdf','ContentType','vector');
end