%TODO brain region filter, olfaction filter.
function [fh,out]=get_olfaction_decoding(sel_meta,opt)
arguments
    sel_meta
    opt.new_data (1,1) logical=true
    %     opt.plot_PCA (1,1) logical=false
    opt.calc_dec (1,1) logical=true
    opt.plot_dec (1,1) logical=true
    opt.ranksum_stats (1,1) logical =false
    opt.lblidx (1,1) double {mustBeMember(opt.lblidx,[5,8])} = 8 %5 for sample 8 for duration
    opt.rpt (1,1) double {mustBeInteger,mustBePositive} = 100
    opt.n_su (1,:) double {mustBeInteger,mustBePositive} =12:12:48;
    opt.type_tag={'type_A','type_B'}
end

if opt.lblidx==5
    behav_tags=[4,8];
    dec_tag="olf";
else
    behav_tags=[3,6];
    dec_tag="dur";
end


%% gen data

if opt.new_data
    decode_mat=struct();
    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
    homedir=ephys.util.getHomedir('type','raw');

%     [dur_sense_mix,dur_indep,~,sens_exclu]=ephys.get_dul_sel('ranksum_stats',opt.ranksum_stats);
    typeAsel=(sel_meta.typeAsel);
    typeBsel=(sel_meta.typeBsel);

    regsel=ismember(meta.reg_tree(1,:),{'CH','BS'}).';

    for ii=reshape(cell2mat(sessmap.keys()),1,[])
        disp(ii)
        fpath=fullfile(homedir,sessmap(ii),"FR_All_1000.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');

        sesssel=meta.sess==ii;

        %         fr_t_align=nan(size(fr,1),size(fr,2)); %trial, su, pre-test bin-averaged
        fr_t_align=mean(fr(:,:,5:7),3); % early/late delay
        %         fr_t_align(trials(:,8)==6,:)=mean(fr(trials(:,8)==6,:,17:28),3); % early/late delay

        decode_mat.(sprintf('s%d',ii)).trials=trials;
        decode_mat.(sprintf('s%d',ii)).fr_type_B=fr_t_align(:,regsel(sesssel) & typeBsel(sesssel));
        decode_mat.(sprintf('s%d',ii)).fr_type_A=fr_t_align(:,regsel(sesssel) & typeAsel(sesssel));
%         decode_mat.(sprintf('s%d',ii)).fr_sens=fr_t_align(:,regsel(sesssel) & sens_exclu(sesssel));
    end

    frmap=struct();
    frmap.(sprintf('fr%d',behav_tags(1))).type_A=containers.Map('KeyType','char','ValueType','any');
    frmap.(sprintf('fr%d',behav_tags(2))).type_A=containers.Map('KeyType','char','ValueType','any');
    frmap.(sprintf('fr%d_e',behav_tags(1))).type_A=containers.Map('KeyType','char','ValueType','any');
    frmap.(sprintf('fr%d_e',behav_tags(2))).type_A=containers.Map('KeyType','char','ValueType','any');
    frmap.(sprintf('fr%d',behav_tags(1))).type_B=containers.Map('KeyType','char','ValueType','any');
    frmap.(sprintf('fr%d',behav_tags(2))).type_B=containers.Map('KeyType','char','ValueType','any');
    frmap.(sprintf('fr%d_e',behav_tags(1))).type_B=containers.Map('KeyType','char','ValueType','any');
    frmap.(sprintf('fr%d_e',behav_tags(2))).type_B=containers.Map('KeyType','char','ValueType','any');

    [dec_c_trlCount,dec_e_trlCount]=deal([]);
    for skey=reshape(fieldnames(decode_mat),1,[])
        trls=decode_mat.(skey{1}).trials;
        csel=trls(:,9)~=0 & trls(:,10)~=0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        esel=trls(:,10)==0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        dec_c_lbl=trls(csel,opt.lblidx);
        dec_e_lbl=trls(esel,opt.lblidx);
        type_A_cfrs=decode_mat.(skey{1}).fr_type_A(csel,:); % [nTrl,nSU]
        type_A_efrs=decode_mat.(skey{1}).fr_type_A(esel,:); % [nTrl,nSU]
        type_B_cfrs=decode_mat.(skey{1}).fr_type_B(csel,:); % [nTrl,nSU]
        type_B_efrs=decode_mat.(skey{1}).fr_type_B(esel,:); % [nTrl,nSU]

        [dur_c_trlSess,dur_e_trlSess]=deal([0,0]);

        for behav_idx=1:2
            dur_c_trlSess(behav_idx)=nnz(dec_c_lbl==behav_tags(behav_idx));
            dur_e_trlSess(behav_idx)=nnz(dec_e_lbl==behav_tags(behav_idx));
            for suidx=1:size(type_A_cfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('fr%d',behav_tags(behav_idx))).type_A(sprintf('%su%d',skey{1},suidx))=type_A_cfrs(dec_c_lbl==behav_tags(behav_idx),suidx);
                frmap.(sprintf('fr%d_e',behav_tags(behav_idx))).type_A(sprintf('%su%d',skey{1},suidx))=type_A_efrs(dec_e_lbl==behav_tags(behav_idx),suidx);
            end
            for suidx=1:size(type_B_cfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('fr%d',behav_tags(behav_idx))).type_B(sprintf('%su%d',skey{1},suidx))=type_B_cfrs(dec_c_lbl==behav_tags(behav_idx),suidx);
                frmap.(sprintf('fr%d_e',behav_tags(behav_idx))).type_B(sprintf('%su%d',skey{1},suidx))=type_B_efrs(dec_e_lbl==behav_tags(behav_idx),suidx);
            end
        end
        dec_c_trlCount=[dec_c_trlCount;dur_c_trlSess];
        dec_e_trlCount=[dec_e_trlCount;dur_e_trlSess];
    end
    save('decode_fr.mat','frmap','dec_c_trlCount','dec_e_trlCount');
else
    load('decode_fr.mat','frmap','dec_c_trlCount','dec_e_trlCount');
end
%% available trials for least repeated condition
%figure();histogram(trlCount(:),1:32)


%% decoding
if opt.calc_dec
    min_trl=20;
    min_e_trl=2;
    su_grp=["type_A","type_B"];
    lbl_grp=["fr"+num2str(behav_tags(1)),"fr"+num2str(behav_tags(2))];
    trl_cnt_grp={dec_c_trlCount,dec_e_trlCount};
    out=struct();

    trl_cnt=trl_cnt_grp;
    lbls=lbl_grp;
    for currgrp=su_grp
        curr_keys=frmap.(char(lbls(1))).(currgrp).keys();
        dur_trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),curr_keys),find(min(trl_cnt{1},[],2)>min_trl & min(trl_cnt{2},[],2)>=min_e_trl));
        curr_keys=curr_keys(dur_trlSel);

        for n_su=opt.n_su
            [result,shuf,result_e]=deal([]);
            for resamp_rpt=1:opt.rpt%15
                sukeys=datasample(curr_keys,n_su,'replace',false);
                rawmat=[...
                    cellfun(@(x) min(x),frmap.(lbls(1)).(currgrp).values(sukeys));...
                    cellfun(@(x) max(x),frmap.(lbls(1)).(currgrp).values(sukeys));...
                    cellfun(@(x) min(x),frmap.(lbls(2)).(currgrp).values(sukeys));...
                    cellfun(@(x) max(x),frmap.(lbls(2)).(currgrp).values(sukeys))];

                [normmat,C,S]=normalize(rawmat,1,"range");
                %              [coeff,~,~]=pca(normmat);

                cv=cvpartition(min_trl,'KFold',10);
                dmat=struct();
                for trlType=lbls
                    dmat.(trlType).(currgrp)=cell2mat(cellfun(@(x) datasample(frmap.(trlType).(currgrp)(x),min_trl),sukeys,'UniformOutput',false));
                    dmat.(strjoin([trlType,"_e"],'')).(currgrp)=cell2mat(cellfun(@(x) datasample(frmap.(strjoin([trlType,"_e"],'')).(currgrp)(x),min_e_trl),sukeys,'UniformOutput',false));
                end

                for kf=1:cv.NumTestSets
                    [xx_train,yy_train,xx_test,yy_test,xx_e_test,yy_e_test]=deal([]);
                    for trlType=lbls
                        xx_train=[xx_train;dmat.(trlType).(currgrp)(training(cv,kf),:)];
                        xx_test=[xx_test;dmat.(trlType).(currgrp)(test(cv,kf),:)];
                        xx_e_test=[xx_e_test;dmat.(strjoin([trlType,"_e"],'')).(currgrp)];
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
            out.(dec_tag).(sprintf('result_%dsu',n_su)).(currgrp)=result;
            out.(dec_tag).(sprintf('shuf_%dsu',n_su)).(currgrp)=shuf;
            out.(dec_tag).(sprintf('etrial_%dsu',n_su)).(currgrp)=result_e;
        end
    end
    save('2type_decoding.mat','out');
elseif opt.plot_dec
    load('2type_decoding.mat','out');
end

if opt.plot_dec
    colors=containers.Map(["result","shuf","etrial"],{'r','k','b'});
    
    su_grp=["type_A","type_B"];
%     lbl_grp={["d3","d6"]};
    lsgrp=["-","--"];
    
    fh=figure('Color','w','Position',[100,100,250,225]);
    hold on;
    for grpidx=1:2
        for ptype=["shuf","etrial","result"]
            %     datamat=cell2mat(arrayfun(@(x) out.(sprintf('%s_%dsu',ptype,x)),[50,100,200,500,1000],'UniformOutput',false));
            n_su=opt.n_su;
            phat=nan(1,numel(n_su));
            pci=nan(2,numel(n_su));
            for nidx=1:numel(n_su)
                dd=out.(dec_tag).(sprintf('%s_%dsu',ptype,n_su(nidx))).(su_grp(grpidx));
                [phat(nidx),pci(:,nidx)]=binofit(nnz(dd),numel(dd));
            end
            fill([n_su,fliplr(n_su)],[pci(1,:),fliplr(pci(2,:))],colors(ptype),'EdgeColor','none','FaceAlpha',0.1);
            ph(grpidx)=plot(n_su,phat,'Color',colors(ptype),'LineStyle',lsgrp(grpidx));
        end
    end
    legend(ph,opt.type_tag,'Location','northoutside','Orientation','horizontal')
    xlabel('Number of neurons')
    ylabel('Classification accuracy (%)')
    ylim([0,1])
    set(gca(),'YTick',0:0.25:1,'YTickLabel',0:25:100)
    xlim([min(n_su),max(n_su)])
    title(dec_tag);
    
%     exportgraphics(fh,'duration_decoding.pdf','ContentType','vector');
end