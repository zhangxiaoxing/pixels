%TODO brain region filter, olfaction filter.
function out=get_duration_decoding(opt)
arguments
    opt.new_data (1,1) logical=true
    %     opt.plot_PCA (1,1) logical=false
    opt.calc_dec (1,1) logical=true
    opt.plot_dec (1,1) logical=true
    opt.ranksum_stats (1,1) logical =false

end
%% gen data

if opt.new_data
    decode_mat=struct();
    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta();
    homedir=ephys.util.getHomedir('type','raw');

    [dur_sense_mix,dur_exclu,~,sens_exclu]=ephys.get_dul_sel('ranksum_stats',opt.ranksum_stats);

    regsel=ismember(meta.reg_tree(1,:),{'CH','BS'}).';

    for ii=reshape(cell2mat(sessmap.keys()),1,[])
        disp(ii)
        fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');

        sesssel=meta.sess==ii;

        %         fr_t_align=nan(size(fr,1),size(fr,2)); %trial, su, pre-test bin-averaged
        fr_t_align=mean(fr(:,:,17:28),3); % early/late delay
        %         fr_t_align(trials(:,8)==6,:)=mean(fr(trials(:,8)==6,:,17:28),3); % early/late delay

        decode_mat.(sprintf('s%d',ii)).trials=trials;
        decode_mat.(sprintf('s%d',ii)).fr_dur=fr_t_align(:,regsel(sesssel) & dur_exclu(sesssel));
        decode_mat.(sprintf('s%d',ii)).fr_sens_dur=fr_t_align(:,regsel(sesssel) & dur_sense_mix(sesssel));
        decode_mat.(sprintf('s%d',ii)).fr_sens=fr_t_align(:,regsel(sesssel) & sens_exclu(sesssel));
    end

    frmap=struct();
    frmap.d3.sensdur=containers.Map('KeyType','char','ValueType','any');
    frmap.d6.sensdur=containers.Map('KeyType','char','ValueType','any');
    frmap.d3_e.sensdur=containers.Map('KeyType','char','ValueType','any');
    frmap.d6_e.sensdur=containers.Map('KeyType','char','ValueType','any');
    frmap.d3.dur=containers.Map('KeyType','char','ValueType','any');
    frmap.d6.dur=containers.Map('KeyType','char','ValueType','any');
    frmap.d3_e.dur=containers.Map('KeyType','char','ValueType','any');
    frmap.d6_e.dur=containers.Map('KeyType','char','ValueType','any');

    frmap.S1.sensdur=containers.Map('KeyType','char','ValueType','any');
    frmap.S2.sensdur=containers.Map('KeyType','char','ValueType','any');
    frmap.S1_e.sensdur=containers.Map('KeyType','char','ValueType','any');
    frmap.S2_e.sensdur=containers.Map('KeyType','char','ValueType','any');

    frmap.S1.sens=containers.Map('KeyType','char','ValueType','any');
    frmap.S2.sens=containers.Map('KeyType','char','ValueType','any');
    frmap.S1_e.sens=containers.Map('KeyType','char','ValueType','any');
    frmap.S2_e.sens=containers.Map('KeyType','char','ValueType','any');

    [dur_c_trlCount,dur_e_trlCount,sens_c_trlCount,sens_e_trlCount]=deal([]);
    for skey=reshape(fieldnames(decode_mat),1,[])
        trls=decode_mat.(skey{1}).trials;
        csel=trls(:,9)~=0 & trls(:,10)~=0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        esel=trls(:,10)==0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        dur_clbl=trls(csel,8);
        dur_elbl=trls(esel,8);
        dur_sens_durcfrs=decode_mat.(skey{1}).fr_sens_dur(csel,:); % [nTrl,nSU]
        dur_sens_durefrs=decode_mat.(skey{1}).fr_sens_dur(esel,:); % [nTrl,nSU]
        dur_durcfrs=decode_mat.(skey{1}).fr_dur(csel,:); % [nTrl,nSU]
        dur_durefrs=decode_mat.(skey{1}).fr_dur(esel,:); % [nTrl,nSU]

        [dur_c_trlSess,dur_e_trlSess,sens_c_trlSess,sens_e_trlSess]=deal([0,0]);
        for trl_samp=[3,6]
            dur_c_trlSess(trl_samp/3)=nnz(dur_clbl==trl_samp);
            dur_e_trlSess(trl_samp/3)=nnz(dur_elbl==trl_samp);
            for suidx=1:size(dur_sens_durcfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('d%d',trl_samp)).sensdur(sprintf('%su%d',skey{1},suidx))=dur_sens_durcfrs(dur_clbl==trl_samp,suidx);
                frmap.(sprintf('d%d_e',trl_samp)).sensdur(sprintf('%su%d',skey{1},suidx))=dur_sens_durefrs(dur_elbl==trl_samp,suidx);
            end
            for suidx=1:size(dur_durcfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('d%d',trl_samp)).dur(sprintf('%su%d',skey{1},suidx))=dur_durcfrs(dur_clbl==trl_samp,suidx);
                frmap.(sprintf('d%d_e',trl_samp)).dur(sprintf('%su%d',skey{1},suidx))=dur_durefrs(dur_elbl==trl_samp,suidx);
            end
        end
        dur_c_trlCount=[dur_c_trlCount;dur_c_trlSess];
        dur_e_trlCount=[dur_e_trlCount;dur_e_trlSess];

        %WIP
        sens_clbl=trls(csel,5);
        sens_elbl=trls(esel,5);
        sens_sens_durcfrs=decode_mat.(skey{1}).fr_sens_dur(csel,:); % [nTrl,nSU]
        sens_sens_durefrs=decode_mat.(skey{1}).fr_sens_dur(esel,:); % [nTrl,nSU]
        sens_senscfrs=decode_mat.(skey{1}).fr_sens(csel,:); % [nTrl,nSU]
        sens_sensefrs=decode_mat.(skey{1}).fr_sens(esel,:); % [nTrl,nSU]

        for trl_samp=1:2
            sens_c_trlSess(trl_samp)=nnz(sens_clbl==trl_samp.*4);
            sens_e_trlSess(trl_samp)=nnz(sens_elbl==trl_samp.*4);
            for suidx=1:size(sens_sens_durcfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('S%d',trl_samp)).sensdur(sprintf('%su%d',skey{1},suidx))=sens_sens_durcfrs(sens_clbl==trl_samp.*4,suidx);
                frmap.(sprintf('S%d_e',trl_samp)).sensdur(sprintf('%su%d',skey{1},suidx))=sens_sens_durefrs(sens_elbl==trl_samp.*4,suidx);
            end
            for suidx=1:size(sens_senscfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('S%d',trl_samp)).sens(sprintf('%su%d',skey{1},suidx))=sens_senscfrs(sens_clbl==trl_samp.*4,suidx);
                frmap.(sprintf('S%d_e',trl_samp)).sens(sprintf('%su%d',skey{1},suidx))=sens_sensefrs(sens_elbl==trl_samp.*4,suidx);
            end
        end
        sens_c_trlCount=[sens_c_trlCount;sens_c_trlSess];
        sens_e_trlCount=[sens_e_trlCount;sens_e_trlSess];

    end
    save('decode_fr.mat','frmap','dur_c_trlCount','dur_e_trlCount','sens_c_trlCount','sens_e_trlCount');
else
    load('decode_fr.mat','frmap','dur_c_trlCount','dur_e_trlCount','sens_c_trlCount','sens_e_trlCount');
end
%% available trials for least repeated condition
%figure();histogram(trlCount(:),1:32)


%% decoding
if opt.calc_dec
    min_trl=20;
    min_e_trl=2;
    su_grp={["sensdur","dur"],["sensdur","sens"]};
    lbl_grp={["d3","d6"],["S1","S2"]};
    dec_tag=["dur","sens"];
    trl_cnt_grp={{dur_c_trlCount,dur_e_trlCount},{sens_c_trlCount,sens_e_trlCount}};
    out=struct();
    for grp_id=1:2
        trl_cnt=trl_cnt_grp{grp_id};
        lbls=lbl_grp{grp_id};
        for currgrp=su_grp{grp_id}
            curr_keys=frmap.(lbls(1)).(currgrp).keys();
            dur_trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),curr_keys),find(min(trl_cnt{1},[],2)>min_trl & min(trl_cnt{2},[],2)>=min_e_trl));
            curr_keys=curr_keys(dur_trlSel);

            for n_su=[5,10,50,100,500]
                [result,shuf,result_e]=deal([]);
                for resamp_rpt=1:50%15
                    sukeys=datasample(curr_keys,n_su);
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
                out.(dec_tag(grp_id)).(sprintf('result_%dsu',n_su)).(currgrp)=result;
                out.(dec_tag(grp_id)).(sprintf('shuf_%dsu',n_su)).(currgrp)=shuf;
                out.(dec_tag(grp_id)).(sprintf('etrial_%dsu',n_su)).(currgrp)=result_e;
            end
        end
    end
    save('duration_decoding.mat','out');
elseif opt.plot_dec
    load('duration_decoding.mat','out');
end

if opt.plot_dec
    colors=containers.Map(["result","shuf","etrial"],{'r','k','b'});
    dec_tag=["dur","sens"];
    su_grp={["sensdur","dur"],["sensdur","sens"]};
    lbl_grp={["d3","d6"],["S1","S2"]};
    lsgrp=["-","--"];
    for dec_idx=1:2
    fh=figure('Color','w','Position',[100,100,250,225]);
    hold on;
    for grpidx=1:2
        for ptype=["shuf","etrial","result"]
            %     datamat=cell2mat(arrayfun(@(x) out.(sprintf('%s_%dsu',ptype,x)),[50,100,200,500,1000],'UniformOutput',false));
            n_su=[5,10,50,100,500];
            phat=nan(1,numel(n_su));
            pci=nan(2,numel(n_su));
            for nidx=1:numel(n_su)
                dd=out.(dec_tag(dec_idx)).(sprintf('%s_%dsu',ptype,n_su(nidx))).(su_grp{dec_idx}(grpidx));
                [phat(nidx),pci(:,nidx)]=binofit(nnz(dd),numel(dd));
            end
            fill([n_su,fliplr(n_su)],[pci(1,:),fliplr(pci(2,:))],colors(ptype),'EdgeColor','none','FaceAlpha',0.1);
            ph(grpidx)=plot(n_su,phat,'Color',colors(ptype),'LineStyle',lsgrp(grpidx));
        end
    end
    legend(ph,{'Mixed modality','Single modality'},'Location','northoutside','Orientation','horizontal')
    xlabel('Number of neurons')
    ylabel('Classification accuracy (%)')
    ylim([0,1])
    set(gca(),'YTick',0:0.1:1,'YTickLabel',0:10:100)
    xlim([0,max(n_su)])
    title(dec_tag(dec_idx));
    end
%     exportgraphics(fh,'duration_decoding.pdf','ContentType','vector');
end