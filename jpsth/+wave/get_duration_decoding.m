%TODO brain region filter, olfaction filter.
function [fh,out]=get_duration_decoding(sel_meta,opt)
arguments
    sel_meta
    opt.new_data (1,1) logical=true
    %     opt.plot_PCA (1,1) logical=false
    opt.calc_dec (1,1) logical=true
    opt.plot_dec (1,1) logical=true
    opt.ranksum_stats (1,1) logical =false
    opt.rpt (1,1) double {mustBeInteger,mustBePositive} = 100

end
%% gen data

if opt.new_data
    decode_mat=struct();
    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta('skip_stats',true);
    homedir=ephys.util.getHomedir();

%     [dur_sense_mix,dur_indep,~,sens_exclu]=ephys.get_dul_sel('ranksum_stats',opt.ranksum_stats);
    dur_dep=ismember(sel_meta.wave_id,1:4);
    dur_indep=ismember(sel_meta.wave_id,5:6);

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
        decode_mat.(sprintf('s%d',ii)).fr_dur_indep=fr_t_align(:,regsel(sesssel) & dur_indep(sesssel));
        decode_mat.(sprintf('s%d',ii)).fr_dur_dep=fr_t_align(:,regsel(sesssel) & dur_dep(sesssel));
%         decode_mat.(sprintf('s%d',ii)).fr_sens=fr_t_align(:,regsel(sesssel) & sens_exclu(sesssel));
    end

    frmap=struct();
    frmap.d3.dur_dep=containers.Map('KeyType','char','ValueType','any');
    frmap.d6.dur_dep=containers.Map('KeyType','char','ValueType','any');
    frmap.d3_e.dur_dep=containers.Map('KeyType','char','ValueType','any');
    frmap.d6_e.dur_dep=containers.Map('KeyType','char','ValueType','any');
    frmap.d3.dur_indep=containers.Map('KeyType','char','ValueType','any');
    frmap.d6.dur_indep=containers.Map('KeyType','char','ValueType','any');
    frmap.d3_e.dur_indep=containers.Map('KeyType','char','ValueType','any');
    frmap.d6_e.dur_indep=containers.Map('KeyType','char','ValueType','any');

    [dur_c_trlCount,dur_e_trlCount]=deal([]);
    for skey=reshape(fieldnames(decode_mat),1,[])
        trls=decode_mat.(skey{1}).trials;
        csel=trls(:,9)~=0 & trls(:,10)~=0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        esel=trls(:,10)==0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        dur_clbl=trls(csel,8);
        dur_elbl=trls(esel,8);
        dur_dep_cfrs=decode_mat.(skey{1}).fr_dur_dep(csel,:); % [nTrl,nSU]
        dur_dep_efrs=decode_mat.(skey{1}).fr_dur_dep(esel,:); % [nTrl,nSU]
        dur_indep_cfrs=decode_mat.(skey{1}).fr_dur_indep(csel,:); % [nTrl,nSU]
        dur_indep_efrs=decode_mat.(skey{1}).fr_dur_indep(esel,:); % [nTrl,nSU]

        [dur_c_trlSess,dur_e_trlSess]=deal([0,0]);
        for trl_samp=[3,6]
            dur_c_trlSess(trl_samp/3)=nnz(dur_clbl==trl_samp);
            dur_e_trlSess(trl_samp/3)=nnz(dur_elbl==trl_samp);
            for suidx=1:size(dur_dep_cfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('d%d',trl_samp)).dur_dep(sprintf('%su%d',skey{1},suidx))=dur_dep_cfrs(dur_clbl==trl_samp,suidx);
                frmap.(sprintf('d%d_e',trl_samp)).dur_dep(sprintf('%su%d',skey{1},suidx))=dur_dep_efrs(dur_elbl==trl_samp,suidx);
            end
            for suidx=1:size(dur_indep_cfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('d%d',trl_samp)).dur_indep(sprintf('%su%d',skey{1},suidx))=dur_indep_cfrs(dur_clbl==trl_samp,suidx);
                frmap.(sprintf('d%d_e',trl_samp)).dur_indep(sprintf('%su%d',skey{1},suidx))=dur_indep_efrs(dur_elbl==trl_samp,suidx);
            end
        end
        dur_c_trlCount=[dur_c_trlCount;dur_c_trlSess];
        dur_e_trlCount=[dur_e_trlCount;dur_e_trlSess];
    end
    save('decode_fr.mat','frmap','dur_c_trlCount','dur_e_trlCount');
else
    load('decode_fr.mat','frmap','dur_c_trlCount','dur_e_trlCount');
end
%% available trials for least repeated condition
%figure();histogram(trlCount(:),1:32)


%% decoding
if opt.calc_dec
    min_trl=20;
    min_e_trl=2;
    su_grp={["dur_dep","dur_indep"]};
    lbl_grp={["d3","d6"]};
    dec_tag=["dur"];
    trl_cnt_grp={{dur_c_trlCount,dur_e_trlCount}};
    out=struct();
    for grp_id=1
        trl_cnt=trl_cnt_grp{grp_id};
        lbls=lbl_grp{grp_id};
        for currgrp=su_grp{grp_id}
            curr_keys=frmap.(lbls(1)).(currgrp).keys();
            dur_trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),curr_keys),find(min(trl_cnt{1},[],2)>min_trl & min(trl_cnt{2},[],2)>=min_e_trl));
            curr_keys=curr_keys(dur_trlSel);

            for n_su=12:12:48
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
    dec_tag=["dur"];
    su_grp={["dur_indep","dur_dep"]};
    lbl_grp={["d3","d6"]};
    lsgrp=["-","--"];
    for dec_idx=1
    fh=figure('Color','w','Position',[100,100,250,225]);
    hold on;
    for grpidx=1:2
        for ptype=["shuf","etrial","result"]
            %     datamat=cell2mat(arrayfun(@(x) out.(sprintf('%s_%dsu',ptype,x)),[50,100,200,500,1000],'UniformOutput',false));
            n_su=12:12:48;
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
    legend(ph,{'Odor independent','Odor dependent'},'Location','northoutside','Orientation','horizontal')
    xlabel('Number of neurons')
    ylabel('Classification accuracy (%)')
    ylim([0,1])
    set(gca(),'YTick',0:0.25:1,'YTickLabel',0:25:100)
    xlim([min(n_su),max(n_su)])
    title(dec_tag(dec_idx));
    end
%     exportgraphics(fh,'duration_decoding.pdf','ContentType','vector');
end