% bin_trial_count=quickStats(ring_list);
midx=2;
rpts=2
if true
inact=false;
cvcorr=cell(1,9);
shufcorr=cell(1,9);
errcorr=cell(1,9);
trial_thres=45;
err_thres=4;
for bin=1:6
    cvcorr{bin+3}=[];
    shufcorr{bin+3}=[];
    errcorr{bin+3}=[];
    for rpt=1:rpts
        disp([bin, rpt]);
        [s1,s2,s1err,s2err]=collect_data(midx,inact,bin,trial_thres,err_thres);
        y=[zeros(trial_thres,1);ones(trial_thres,1)];
        cv=cvpartition(trial_thres,'KFold',10);
        for kf=1:cv.NumTestSets
            s1kf=s1(training(cv,kf),:);
            s2kf=s2(training(cv,kf),:);
            varsel=var(s1kf,0,1)>0 & var(s2kf,0,1)>0;
            Xkf=[s1kf(:,varsel);s2kf(:,varsel)];
            Xerr=[s1err(:,varsel);s2err(:,varsel)];
            yerr=[zeros(size(s1err,1),1);ones(size(s2err,1),1)];
            ykf=y([training(cv,kf);training(cv,kf)]);
            CVLDAModel=fitcdiscr(Xkf,ykf,'DiscrimType','linear');
            s1Tkf=s1(test(cv,kf),:);
            s2Tkf=s2(test(cv,kf),:);
            XTkf=[s1Tkf(:,varsel);s2Tkf(:,varsel)];
            yTkf=y([test(cv,kf);test(cv,kf)]);
            yshufTkf=yTkf(randperm(numel(yTkf)));
            cvresult=CVLDAModel.predict(XTkf)==yTkf;
            cvshufresult=CVLDAModel.predict(XTkf)==yshufTkf;
            cv_err_result=CVLDAModel.predict(Xerr)==yerr;
            cvcorr{bin+3}=[cvcorr{bin+3};cvresult];
            shufcorr{bin+3}=[shufcorr{bin+3};cvshufresult];
            errcorr{bin+3}=[errcorr{bin+3};cv_err_result];
        end
    end
end
fh=figure('Color','w','Position',[100,100,195,160]);
hold on
lossci=(cell2mat(cellfun(@(x) bootci(1000,@(y) mean(y), x),cvcorr([4:9]),'UniformOutput',false)))*100;
shufci=(cell2mat(cellfun(@(x) bootci(1000,@(y) mean(y), x),shufcorr([4:9]),'UniformOutput',false)))*100;
errci=(cell2mat(cellfun(@(x) bootci(1000,@(y) mean(y), x),errcorr([4:9]),'UniformOutput',false)))*100;

fill([1:6,fliplr(1:6)],[lossci(1,:),fliplr(lossci(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
fill([1:6,fliplr(1:6)],[shufci(1,:),fliplr(shufci(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
fill([1:6,fliplr(1:6)],[errci(1,:),fliplr(errci(2,:))],'b','EdgeColor','none','FaceAlpha',0.2)

plot([1:6],(cellfun(@(x) mean(x),cvcorr([4:9])))*100,'-r');
plot([1:6],(cellfun(@(x) mean(x),shufcorr([4:9])))*100,'-k');
plot([1:6],(cellfun(@(x) mean(x),errcorr([4:9])))*100,'-k');
xlabel('Time (s)')
ylabel('Sample classification accuracy (%)')
xlim([-2,6])
set(gca,'XTick',-1:2:5,'XTickLabel',{'ITI','1','3','5'})
% exportgraphics(fh,'congru_3ring_decode.pdf','ContentType','vector')
savefig(fh,sprintf('congru_%dring_decode.fig',midx+2),'compact');
print(fh,sprintf('congru_%dring_decode.png',midx+2),'-dpng','-r300');

end
act_cvcorr=cvcorr;
act_cvshuf=shufcorr;

if false
inact=true;
cvcorr=cell(1,9);
shufcorr=cell(1,9);
for bin=[-2,1:6]
    cvcorr{bin+3}=[];
    shufcorr{bin+3}=[];
    for rpt=1:rpts
        disp([bin, rpt]);
        [s1,s2]=collect_data(midx,inact,bin,trial_thres,err_thres);
        y=[zeros(trial_thres,1);ones(trial_thres,1)];
        cv=cvpartition(trial_thres,'KFold',10);
        for kf=1:cv.NumTestSets
            s1kf=s1(training(cv,kf),:);
            s2kf=s2(training(cv,kf),:);
            varsel=var(s1kf,0,1)>0 & var(s2kf,0,1)>0;
            Xkf=[s1kf(:,varsel);s2kf(:,varsel)];
            ykf=y([training(cv,kf);training(cv,kf)]);
            CVLDAModel=fitcdiscr(Xkf,ykf,'DiscrimType','linear');
            s1Tkf=s1(test(cv,kf),:);
            s2Tkf=s2(test(cv,kf),:);
            XTkf=[s1Tkf(:,varsel);s2Tkf(:,varsel)];
            yTkf=y([test(cv,kf);test(cv,kf)]);
            yshufTkf=yTkf(randperm(numel(yTkf)));
            cvresult=CVLDAModel.predict(XTkf)==yTkf;
            cvshufresult=CVLDAModel.predict(XTkf)==yshufTkf;
            cvcorr{bin+3}=[cvcorr{bin+3};cvresult];
            shufcorr{bin+3}=[shufcorr{bin+3};cvshufresult];
        end
    end
end
inact_cvcorr=cvcorr;
inact_cvshuf=shufcorr;


save('rings_decode.mat','act_cvcorr','act_cvshuf','inact_cvcorr','inact_cvshuf')

fh=figure('Color','w','Position',[100,100,195,160]);
hold on
lossci=(cell2mat(cellfun(@(x) bootci(1000,@(y) mean(y), x),inact_cvcorr([1,4:9]),'UniformOutput',false)))*100;
shufci=(cell2mat(cellfun(@(x) bootci(1000,@(y) mean(y), x),inact_cvshuf([1,4:9]),'UniformOutput',false)))*100;

fill([-1,1:6,fliplr(1:6),-1],[lossci(1,:),fliplr(lossci(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
fill([-1,1:6,fliplr(1:6),-1],[shufci(1,:),fliplr(shufci(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)

plot([-1,1:6],(cellfun(@(x) mean(x),inact_cvcorr([1,4:9])))*100,'-r');
plot([-1,1:6],(cellfun(@(x) mean(x),inact_cvshuf([1,4:9])))*100,'-k');
xlabel('Time (s)')
ylabel('Sample classification accuracy (%)')
xlim([-2,6])
set(gca,'XTick',-1:2:5,'XTickLabel',{'ITI','1','3','5'})
% exportgraphics(fh,'congru_3ring_inact_decode.pdf','ContentType','vector')

savefig(fh,sprintf('congru_%dring_inact_decode.fig',midx+2),'compact');
print(fh,sprintf('congru_%dring_inact_decode.png',midx+2),'-dpng','-r300');
end
if isunix
    quit(0) 
end


% keyboard()
% 

function [S1_rings_fr,S2_rings_fr,S1_rings_err_fr,S2_rings_err_fr]=collect_data(midx,inact,bin,trials,err_trials)
datapath='ring_freq';
if exist('inact','var') && inact
    flist=dir(fullfile(datapath,sprintf('*freq_inact_%d_*.mat',midx+2)));
else
    flist=dir(fullfile(datapath,sprintf('*freq_%d_*.mat',midx+2)));
end

ring_list=cell(0,8);
err_ring_list=cell(0,8);
for i=1:size(flist,1)
    fn=fullfile(datapath,flist(i).name);
    load(fn,'all_rings','err_rings');
    ring_list=[ring_list;all_rings];
    err_ring_list=[err_ring_list;err_rings];
end

ringOccur=sum(cell2mat(ring_list(:,5:6)),2);
totalTrial=sum(cellfun(@(x) size(x,1), ring_list(:,7:8)),2);

ring_list=ring_list(ringOccur>=totalTrial,:);
err_ring_list=err_ring_list(ringOccur>=totalTrial,:);

S1_rings_fr=[];
S2_rings_fr=[];
S1_rings_err_fr=[];
S2_rings_err_fr=[];
for i=1:size(ring_list,1)
    if size(ring_list{i,7},1)<=trials ...
            || size(ring_list{i,8},1)<=trials ...
            || size(err_ring_list{i,7},1)<=err_trials ...
            || size(err_ring_list{i,8},1)<=err_trials
        continue
    end
    s1=(ring_list{i,7}(:,bin+3));
    s2=(ring_list{i,8}(:,bin+3));
    s1_err=(err_ring_list{i,7}(:,bin+3));
    s2_err=(err_ring_list{i,8}(:,bin+3));
    
    if nnz(~isnan(s1))>trials && nnz(~isnan(s2))>trials...
        && nnz(~isnan(s1_err))>err_trials && nnz(~isnan(s2_err))>err_trials
        s1samp=randsample(s1(~isnan(s1)),trials);
        s2samp=randsample(s2(~isnan(s2)),trials);
        s1samp_err=randsample(s1_err(~isnan(s1_err)),err_trials);
        s2samp_err=randsample(s2_err(~isnan(s2_err)),err_trials);
        S1_rings_fr=[S1_rings_fr,s1samp];
        S2_rings_fr=[S2_rings_fr,s2samp];
        S1_rings_err_fr=[S1_rings_err_fr,s1samp_err];
        S2_rings_err_fr=[S2_rings_err_fr,s2samp_err];
    end
end
end




function out=quickStats(ring_list)
out=[];
for bin=1:6
    for i=1:size(ring_list,1)
        s1t=nnz(~isnan(ring_list{i,7}(:,bin+3)));
        s2t=nnz(~isnan(ring_list{i,8}(:,bin+3)));
        if s1t>0 && s2t>0
            out=[out;bin,i,s1t,s2t];
        end
    end
end

for bin=1:6
    disp(nnz(out(:,1)==bin & out(:,3)>30 & out(:,4)>30))
end

end
