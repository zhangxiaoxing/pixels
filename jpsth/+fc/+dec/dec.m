function out=dec(decdata,opt)
arguments
    decdata struct
    opt.decoder (1,:) char {mustBeMember(opt.decoder,{'SVM','LDA','NB'})} = 'SVM'
    opt.svmcost (1,1) double = 1
    opt.varsel (1,1) logical = false
    opt.foldN (1,1) double = 2
end

s1=decdata.s1;
s2=decdata.s2;

e1=decdata.e1;
e2=decdata.e2;

bins=size(s1,3);
trlN=size(s1,2);
trlE=size(e1,2);
rpts=size(s1,4);
cv=cvpartition(trlN,'KFold',opt.foldN);
y=[zeros(trlN,1);ones(trlN,1)];

out=struct();
out.cvcorr=cell(1,bins);
out.shufcorr=cell(1,bins);
out.errcorr=cell(1,bins);

%     errcorr{bin+3}=[];

for rpt=1:rpts
    fprintf('rpts %d of %d\n',rpt,rpts);
    for bin=1:bins
        for kf=1:cv.NumTestSets
            s1kf=s1(:,training(cv,kf),bin,rpt);
            s2kf=s2(:,training(cv,kf),bin,rpt);
            Xkf=cat(2,s1kf,s2kf)';
            ykf=y([training(cv,kf);training(cv,kf)]);
            if opt.varsel
                s1Tkf=s1(varsel,test(cv,kf),bin,rpt);
                s2Tkf=s2(varsel,test(cv,kf),bin,rpt);
            else
                s1Tkf=s1(:,test(cv,kf),bin,rpt);
                s2Tkf=s2(:,test(cv,kf),bin,rpt);
            end
            XTkf=cat(2,s1Tkf,s2Tkf)';
            yTkf=y([test(cv,kf);test(cv,kf)]);
            yshufTkf=yTkf(randperm(numel(yTkf)));


            if strcmp(opt.decoder,'SVM')
                SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',false,'Cost',[0,opt.svmcost;opt.svmcost,0]);
                modelPredict=SVMM.predict(XTkf);
            elseif strcmp(opt.decoder,'LDA')
                LDAM=fitcdiscr(Xkf,ykf,'DiscrimType','pseudolinear');
                modelPredict=LDAM.predict(XTkf);
            elseif strcmp(opt.decoder,'NB')
                NBM=fitcnb(Xkf,ykf);
                modelPredict=NBM.predict(XTkf);
            end

            cvresult=modelPredict==yTkf;
            cvshufresult=modelPredict==yshufTkf;

            out.cvcorr{bin}=cat(1,out.cvcorr{bin},cvresult);
            out.shufcorr{bin}=cat(1,out.shufcorr{bin},cvshufresult);
            if strcmp(opt.decoder,'SVM')
                cv_err_result=SVMM.predict([e1(:,:,bin,rpt).';e2(:,:,bin,rpt).'])==((1:2*trlE)>trlE).';
            elseif strcmp(opt.decoder,'LDA')
                cv_err_result=LDAM.predict([e1(:,:,bin,rpt).';e2(:,:,bin,rpt).'])==((1:2*trlE)>trlE).';
            end
            out.errcorr{bin}=cat(1,out.errcorr{bin},cv_err_result);
        end
    end
end
end

function plotSV
sv=SVMM.SupportVectors;
for ii=1:2:200
    fh=figure()
    gscatter(Xkf(:,ii),Xkf(:,ii+1),ykf)
    hold on
    plot(sv(:,ii),sv(:,ii+1),'ko','MarkerSize',10)
    waitfor(fh)
end
end