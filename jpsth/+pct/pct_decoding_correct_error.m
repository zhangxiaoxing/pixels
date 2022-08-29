%TODO brain region filter, olfaction filter.
function out=pct_decoding_correct_error(pct_meta,ids,opt)
arguments
    pct_meta
    ids (1,:) double {mustBeInteger,mustBePositive}
    opt.new_data (1,1) logical=true
    opt.calc_dec (1,1) logical=true
    opt.lblidx (1,1) double {mustBeMember(opt.lblidx,[5,8])} = 8 %5 for sample 8 for duration
    opt.rpt (1,1) double {mustBeInteger,mustBePositive} = 25 % x10 fold x4 test/fold = 1000
    opt.n_su (1,:) double {mustBeInteger,mustBePositive} =100;
end

if opt.lblidx==5
    behav_tags=[4,8];
    dec_tag="olf";
else
    behav_tags=[3,6];
    dec_tag="dur";
end

% select su by wave_id>>>>>>>>>>>>>>>
susel=ismember(pct_meta.wave_id,ids);
%% gen data
if opt.new_data
    decode_mat=struct();
    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
    homedir=ephys.util.getHomedir('type','raw');
    regsel=ismember(meta.reg_tree(1,:),{'CH','BS'}).';
    frmap=struct();
    
    for ii=reshape(cell2mat(sessmap.keys()),1,[])
        disp(ii)
        fpath=fullfile(homedir,sessmap(ii),"FR_All_1000.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        sesssel=meta.sess==ii;
        fr_t_align=mean(fr(:,:,5:7),3);

        decode_mat.(sprintf('s%d',ii)).trials=trials;
        decode_mat.(sprintf('s%d',ii)).per_su=fr_t_align(:,regsel(sesssel) & susel(sesssel));
    end
    frmap.(sprintf('fr%d',behav_tags(1))).per_su=containers.Map('KeyType','char','ValueType','any');
    frmap.(sprintf('fr%d',behav_tags(2))).per_su=containers.Map('KeyType','char','ValueType','any');

    trlCount=[];
    for skey=reshape(fieldnames(decode_mat),1,[])
        trls=decode_mat.(skey{1}).trials;
        csel=trls(:,9)~=0 & trls(:,10)~=0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        esel=trls(:,10)==0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);

        dec_c_lbl=trls(csel,opt.lblidx);
        dec_e_lbl=trls(esel,opt.lblidx);        

        c_e_trlSess=[0,0,0,0];
        for behav_idx=1:2
            c_e_trlSess(behav_idx)=nnz(dec_c_lbl==behav_tags(behav_idx));
            c_e_trlSess(behav_idx+2)=nnz(dec_e_lbl==behav_tags(behav_idx));
            cfrs=decode_mat.(skey{1}).per_su(csel,:); % [nTrl,nSU]
            efrs=decode_mat.(skey{1}).per_su(esel,:); % [nTrl,nSU]
            for suidx=1:size(cfrs,2) % TODO: vectorized batch init
                frmap.(sprintf('fr%d',behav_tags(behav_idx))).per_su(sprintf('c%su%d',skey{1},suidx))=cfrs(dec_c_lbl==behav_tags(behav_idx),suidx);
                frmap.(sprintf('fr%d',behav_tags(behav_idx))).per_su(sprintf('e%su%d',skey{1},suidx))=efrs(dec_e_lbl==behav_tags(behav_idx),suidx);
            end
        end
        trlCount=[trlCount;c_e_trlSess];
    end
    save(sprintf('pct_decode_corr_err_fr_%s.mat',dec_tag),'frmap','trlCount');
else
    load(sprintf('pct_decode_corr_err_fr_%s.mat',dec_tag),'frmap','trlCount');
end
%% available trials for least repeated condition
%figure();histogram(trlCount(:),1:32)


%% decoding
if opt.calc_dec
    min_c_trl=20;
    min_e_trl=2;
    lbls=["fr"+num2str(behav_tags(1)),"fr"+num2str(behav_tags(2))];
    out=struct();

    %%  bidx=1:(numel(ebound)-1) >>>>>>>>>>>>>>>>>>

    curr_keys=frmap.(lbls(1)).per_su.keys();
    c_keys=curr_keys(startsWith(curr_keys,'c'));
    %         e_keys=curr_keys(startsWith(curr_keys,'e'));

    ntrlsel=find(min(trlCount(:,1:2),[],2)>min_c_trl...
        & min(trlCount(:,3:4),[],2)>min_e_trl);

    trl_sel=ismember(cellfun(@(x) ...
        str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),c_keys),ntrlsel);
    c_keys=c_keys(trl_sel);

    for n_su=opt.n_su
        [c_result,e_result]=deal([]);
        for resamp_rpt=1:opt.rpt%15
            c_sukeys=datasample(c_keys,n_su,'replace',false);
            e_sukeys=replace(c_sukeys,'cs','es');
            rawmat=[... % for scaling >>>>>>>>>>
                cellfun(@(x) min(x),frmap.(lbls(1)).per_su.values(c_sukeys));...
                cellfun(@(x) max(x),frmap.(lbls(1)).per_su.values(c_sukeys));...
                cellfun(@(x) min(x),frmap.(lbls(2)).per_su.values(c_sukeys));...
                cellfun(@(x) max(x),frmap.(lbls(2)).per_su.values(c_sukeys))];

            [normmat,C,S]=normalize(rawmat,1,"range");
            % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            %              [coeff,~,~]=pca(normmat);

            cv=cvpartition(min_c_trl,'KFold',10);
            dmat=struct();
            for trlType=lbls
                dmat.(trlType).c=cell2mat(...
                    cellfun(@(x) datasample(frmap.(trlType).per_su(x),min_c_trl),...
                    c_sukeys,'UniformOutput',false));

                dmat.(trlType).e=cell2mat(...
                    cellfun(@(x) datasample(frmap.(trlType).per_su(x),min_e_trl),...
                    e_sukeys,'UniformOutput',false));
            end

            for kf=1:cv.NumTestSets
                [xx_train,yy_train,xx_test,yy_test,xx_err,yy_err]=deal([]);
                for trlType=lbls
                    xx_train=[xx_train;dmat.(trlType).c(training(cv,kf),:)];
                    xx_test=[xx_test;dmat.(trlType).c(test(cv,kf),:)];
                    xx_err=[xx_err;dmat.(trlType).e];
                    yy_train=[yy_train;repmat(trlType,nnz(training(cv,kf)),1)]; %labels
                    yy_test=[yy_test;repmat(trlType,nnz(test(cv,kf)),1)];
                    yy_err=[yy_err;repmat(trlType,size(dmat.(trlType).e,1),1)];
                end
                norm_train=normalize(xx_train,'center',C,'scale',S);
                SVMM=fitcecoc(norm_train,yy_train);
                norm_test=normalize(xx_test,'center',C,'scale',S);
                modelPredict=SVMM.predict(norm_test);
                norm_err=normalize(xx_err,'center',C,'scale',S);
                err_pred=SVMM.predict(norm_err);
                c_result=[c_result;string(modelPredict)==yy_test];
                e_result=[e_result;string(err_pred)==yy_err];
            end
        end
        out.(dec_tag).(sprintf('c_result_%dsu',n_su))=c_result;
        out.(dec_tag).(sprintf('e_result_%dsu',n_su))=e_result;
    end
    %     end<<<<<<<<<<<<<<<<<<<<<<<<<<<
    save(sprintf('correct_err_decoding_%s.mat',dec_tag),'out');
elseif opt.plot_dec
    load(sprintf('correct_err_band_decoding_%s.mat',dec_tag),'out');
end
