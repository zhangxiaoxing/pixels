function svm_meta=comb_coding(opt)
arguments
    opt.hist_model (1,1) logical = true % modeling history effect
    opt.lo (1,1) double {mustBeInteger,mustBePositive} = 1
    opt.hi (1,1) double {mustBeInteger,mustBePositive} = 33028
    opt.save_file (1,1) logical = false
end

su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
if opt.hi>numel(su_meta.allcid)
    opt.hi=numel(su_meta.allcid);
end

svm_meta=struct();
svm_meta.polf=ones(numel(su_meta.allcid),3);
svm_meta.pdur=ones(numel(su_meta.allcid),3);
svm_meta.lo=opt.lo;
svm_meta.hi=opt.hi;
currstem=[];
homedir=ephys.util.getHomedir('type','raw');

for ii=opt.lo:opt.hi
    if ~strcmp(currstem,su_meta.allpath{ii})
        currstem=su_meta.allpath{ii};
        fpath=replace(fullfile(homedir,su_meta.allpath{ii},'FR_All_1000.hdf5'),'\',filesep());
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        fr=h5read(fpath,'/FR_All');
        disp(num2str(ii)+" total trials "+num2str(nnz(trials(:,9) & trials(:,10))))
        s1d3t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==3);
        s2d3t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==3);
        s1d6t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6);
        s2d6t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6);

    end

    oneuid=suid(suid==su_meta.allcid(ii));
    
    if opt.hist_model
        nblk=min([numel(s1d3t),numel(s2d3t),numel(s1d6t),numel(s2d6t)]);
        %% sample
        [decdata.s1,decdata.s2]=deal([]);
        for bb=1:nblk
            decdata.s1=cat(2,decdata.s1,[squeeze(fr(s1d3t(bb),suid==oneuid,5:7))-squeeze(fr(s1d6t(bb),suid==oneuid,5:7))]);
            decdata.s2=cat(2,decdata.s2,[squeeze(fr(s2d3t(bb),suid==oneuid,5:7))-squeeze(fr(s2d6t(bb),suid==oneuid,5:7))]);
        end
        decout=dec(decdata);
%         [~,~,svm_meta.polf(ii)]=crosstab(1:(2*numel(decout.cvcorr))>numel(decout.cvcorr),[decout.cvcorr;decout.shufcorr]);
        svm_meta.polf(ii,:)=[sum(decout.cvcorr),numel(decout.cvcorr),2*(1-binocdf(sum(decout.cvcorr),numel(decout.cvcorr),0.5))];

        %% duration
        [decdata.s1,decdata.s2]=deal([]);
        for bb=1:nblk
            decdata.s1=cat(2,decdata.s1,[squeeze(fr(s1d3t(bb),suid==oneuid,5:7))-squeeze(fr(s2d3t(bb),suid==oneuid,5:7))]);
            decdata.s2=cat(2,decdata.s2,[squeeze(fr(s1d6t(bb),suid==oneuid,5:7))-squeeze(fr(s2d6t(bb),suid==oneuid,5:7))]);
        end
        decout=dec(decdata);
        %         [~,~,svm_meta.pdur(ii)]=crosstab(1:(2*numel(decout.cvcorr))>numel(decout.cvcorr),[decout.cvcorr;decout.shufcorr]);
        svm_meta.pdur(ii,:)=[sum(decout.cvcorr),numel(decout.cvcorr),2*(1-binocdf(sum(decout.cvcorr),numel(decout.cvcorr),0.5))];
    else
%% sample
        decdata.s1=squeeze(fr([s1d3t;s1d6t],suid==oneuid,5:7)).';
        decdata.s2=squeeze(fr([s2d3t;s2d6t],suid==oneuid,5:7)).';
        decout=dec(decdata);
        svm_meta.polf(ii,:)=[sum(decout.cvcorr),numel(decout.cvcorr),2*(1-binocdf(sum(decout.cvcorr),numel(decout.cvcorr),0.5))];
%% duration
        decdata.s1=squeeze(fr([s1d3t;s2d3t],suid==oneuid,5:7)).';
        decdata.s2=squeeze(fr([s1d6t;s2d6t],suid==oneuid,5:7)).';
        decout=dec(decdata);
        svm_meta.pdur(ii,:)=[sum(decout.cvcorr),numel(decout.cvcorr),2*(1-binocdf(sum(decout.cvcorr),numel(decout.cvcorr),0.5))];
    end
    if opt.save_file
        if opt.hist_model, suffix='hist';else, suffix='';end
        fn=sprintf('comb_coding_%05d_%05d_%s.mat',opt.lo,opt.hi,suffix);
        save(fn,'svm_meta');
    end
    

end
end

%%
function out=dec(decdata)
s1=decdata.s1;
s2=decdata.s2;

trlN=min(size(s1,2),size(s2,2));
cv=cvpartition(trlN,'KFold',10);
y=[zeros(trlN,1);ones(trlN,1)];

out=struct();
out.cvcorr=[];
% out.shufcorr=[];

for kf=1:cv.NumTestSets
    s1kf=s1(:,training(cv,kf));
    s2kf=s2(:,training(cv,kf));
    Xkf=cat(2,s1kf,s2kf)';
    ykf=y([training(cv,kf);training(cv,kf)]);
    s1Tkf=s1(:,test(cv,kf));
    s2Tkf=s2(:,test(cv,kf));

    XTkf=cat(2,s1Tkf,s2Tkf)';
    yTkf=y([test(cv,kf);test(cv,kf)]);
%     yshufTkf=yTkf(randperm(numel(yTkf)));

    SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear');
    modelPredict=SVMM.predict(XTkf);

    cvresult=modelPredict==yTkf;
%     cvshufresult=modelPredict==yshufTkf;

    out.cvcorr=cat(1,out.cvcorr,cvresult);
%     out.shufcorr=cat(1,out.shufcorr,cvshufresult);
end
end

