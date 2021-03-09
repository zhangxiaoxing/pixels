function out=dec(decdata,opt)
arguments
    decdata struct
    opt.decoder (1,:) char {mustBeMember(opt.decoder,{'SVM','LDA','NB'})} = 'SVM'
end

s1=decdata.s1;
s2=decdata.s2;

bins=size(s1,3);
trlN=size(s1,2);
rpts=size(s1,4);
cv=cvpartition(trlN,'KFold',10);
y=[zeros(trlN,1);ones(trlN,1)];

out=struct();
out.cvcorr=cell(1,bins);
out.shufcorr=cell(1,bins);

%     errcorr{bin+3}=[];

for rpt=1:rpts
    fprintf('rpts %d of %d\n',rpt,rpts);
    for bin=1:bins
        for kf=1:cv.NumTestSets
            s1kf=s1(:,training(cv,kf),bin,rpt);
            s2kf=s2(:,training(cv,kf),bin,rpt);
            varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
            s1kf=s1kf(varsel,:);
            s2kf=s2kf(varsel,:);
%             Xerr=[s1err(:,varsel);s2err(:,varsel)];
%             yerr=[zeros(size(s1err,1),1);ones(size(s2err,1),1)];
            Xkf=cat(2,s1kf,s2kf)';
            ykf=y([training(cv,kf);training(cv,kf)]);
            s1Tkf=s1(varsel,test(cv,kf),bin,rpt);
            s2Tkf=s2(varsel,test(cv,kf),bin,rpt);
            XTkf=cat(2,s1Tkf,s2Tkf)';
            yTkf=y([test(cv,kf);test(cv,kf)]);
            yshufTkf=yTkf(randperm(numel(yTkf)));

            
            if strcmp(opt.decoder,'SVM')
                SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',true);
                modelPredict=SVMM.predict(XTkf);
            elseif strcmp(opt.decoder,'LDA')
                LDAM=fitcdiscr(Xkf,ykf);
                modelPredict=LDAM.predict(XTkf);
            elseif strcmp(opt.decoder,'NB')
                NBM=fitcnb(Xkf,ykf);
                modelPredict=NBM.predict(XTkf);                
            end
            
            cvresult=modelPredict==yTkf;
            cvshufresult=modelPredict==yshufTkf;
%             cv_err_result=CVLDAModel.predict(Xerr)==yerr;
            out.cvcorr{bin}=cat(1,out.cvcorr{bin},cvresult);
            out.shufcorr{bin}=cat(1,out.shufcorr{bin},cvshufresult);
%             errcorr{bin+3}=[errcorr{bin+3};cv_err_result];
        end
    end
end
end



