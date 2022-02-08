function out=decode_one(datain,opt)
arguments
    datain (:,1) cell
    opt.keyP (1,:) char = 'pref'
    opt.keyN (1,:) char = 'nonp'
    opt.unit_count (1,:) double = [10,50,100,500];
    opt.trial_count (1,1) double = 10
    opt.rpt (1,1) double = 50
    opt.pre_pca (1,1) logical = false
%     opt.datatype (1,:) char {mustBeMember(opt.datatype,{'SU','FCSP','loops'})} = 'SU';
end

min_trl_cnt=cellfun(@(x) min([numel(x(opt.keyP)),numel(x(opt.keyN))]), datain);
%TODO possible error trials criteria
data=datain(min_trl_cnt>opt.trial_count);
datatype='SU';
if isKey(data{1},'meta') && numel(data{1}('meta'))>4
    members=data{1}('meta');
    if numel(members{5})>2
        datatype='loops';
    elseif numel(members{5})==2
        datatype='FCSP';
    end
end

out=cell(0,1);
for n_unit=opt.unit_count
    [result,shuf,result_e,su_count]=deal([]);
    for resamp_rpt=1:opt.rpt
        unit_idces=randsample(size(data,1),n_unit);
        switch datatype
            case 'loops'
                metas=arrayfun(@(x) data{x}('meta'),unit_idces,'UniformOutput',false);
                curr_su_count=numel(unique(cell2mat(cellfun(@(x) x{4}*100000+x{5},metas.','UniformOutput',false))));
            case 'SU'
                curr_su_count=n_unit;
        end

        Pmat=cell2mat(arrayfun(@(x) datasample(data{x}(opt.keyP),opt.trial_count),unit_idces,'UniformOutput',false).');
        Nmat=cell2mat(arrayfun(@(x) datasample(data{x}(opt.keyN),opt.trial_count),unit_idces,'UniformOutput',false).');

        joinmat=[Pmat;Nmat];

        [normmat,C,S]=normalize(joinmat,1,"range");
        if opt.pre_pca
            [coeff,~,~]=pca(normmat);
        end
        cv=cvpartition(opt.trial_count,'KFold',10);

        Pdata=normmat(1:size(normmat,1)/2,:);
        Ndata=normmat(size(normmat,1)/2+1:end,:);
        
        %             dmat.(strjoin([trlType,"_e"],''))=cell2mat(cellfun(@(x) datasample(frmap.(strjoin([trlType,"_e"],''))(x),min_e_trl),unit_ids,'UniformOutput',false));
        
        for kf=1:cv.NumTestSets
            xx_train=[Pdata(training(cv,kf),:);Ndata(training(cv,kf),:)];
            xx_test=[Pdata(test(cv,kf),:);Ndata(test(cv,kf),:)];
%             xx_e_test=[xx_e_test;dmat.(strjoin([trlType,"_e"],''))];
            yy_train=[ones(nnz(training(cv,kf)),1);-1*ones(nnz(training(cv,kf)),1)];
            yy_test=[ones(nnz(test(cv,kf)),1);-1*ones(nnz(test(cv,kf)),1)];
%             yy_e_test=[yy_e_test;repmat(trlType,min_e_trl,1)];

            norm_train=normalize(xx_train,'center',C,'scale',S);
            norm_test=normalize(xx_test,'center',C,'scale',S);
            if opt.pre_pca
                pc_train=norm_train*coeff;
                SVMM=fitcsvm(pc_train,yy_train,'KernelFunction','rbf');
                pc_test=norm_test*coeff;
                modelPredict=SVMM.predict(pc_test);
            else
                SVMM=fitcsvm(norm_train,yy_train,'KernelFunction','rbf');
                modelPredict=SVMM.predict(norm_test);
            end
            result=[result;modelPredict==yy_test];
            shuf=[shuf;randsample(modelPredict,numel(modelPredict))==yy_test];
            su_count=[su_count;curr_su_count;curr_su_count];
%             norm_e_test=normalize(xx_e_test,'center',C,'scale',S);
%             pc_e=norm_e_test*coeff;
%             e_predict=SVMM.predict(pc_e);
%             result_e=[result_e;e_predict==yy_e_test];
        end
    end
    onestep=struct();
    onestep.result=result;
    onestep.shuf=shuf;
    onestep.su_count=su_count;
    onestep.unit_count=n_unit;
    out{end+1,1}=onestep;
end




