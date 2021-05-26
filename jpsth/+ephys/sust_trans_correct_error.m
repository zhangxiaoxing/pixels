function sust_trans_correct_error(type,opt)
arguments
    type (1,:) char {mustBeMember(type,{'sust','trans'})}
    opt.lastbin (1,1) logical = false
end
meta=ephys.util.load_meta();
stats=struct();
figure('Color','w','Position',[100,100,400,400]);
hold on;
fidx=1;
colors={'r','b'};
if strcmp(type,'sust')
    types=[1,3];
else
    types=[2,4];
end
for onetype=types
    typesel=find(meta.mem_type==onetype);
    if strcmp(type,'trans')
        if opt.lastbin
            typesel=typesel(1:5:numel(typesel));
        else
            typesel=typesel(1:50:numel(typesel));
        end
    end
    homedir=ephys.util.getHomedir('type','raw');
    idx=1;
    stats.(sprintf('type%d',onetype))=nan(numel(typesel),4);
    for ii=typesel
        fpath=fullfile(homedir,meta.allpath{ii},'FR_All_1000.hdf5');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        fr=h5read(fpath,'/FR_All');
        if opt.lastbin
            bins=6;
            if meta.per_bin(6,ii)==0
                continue
            end
        else
            bins=find(meta.per_bin(:,ii)~=0)+4;
        end
        
        cs1=mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6,suid==meta.allcid(ii),bins),3);
        cs2=mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6,suid==meta.allcid(ii),bins),3);
        es1=mean(fr(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==6,suid==meta.allcid(ii),bins),3);
        es2=mean(fr(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==6,suid==meta.allcid(ii),bins),3);
        
        delaymm=mean([mean(cs1),mean(cs2)]);
        delaystd=std([cs1;cs2]);
        
        cs1=(cs1-delaymm)./delaystd;
        cs2=(cs2-delaymm)./delaystd;
        es1=(es1-delaymm)./delaystd;
        es2=(es2-delaymm)./delaystd;
        
        ch(fidx)=scatter(mean(cs1),mean(cs2),36,colors{fidx},'MarkerFaceColor',colors{fidx},'MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');
        eh=scatter(mean(es1),mean(es2),16,'k','MarkerFaceColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
        
        plot([mean(cs1),mean(es1)],[mean(cs2),mean(es2)],':','LineWidth',0.5,'Color',colors{fidx});
        stats.(sprintf('type%d',onetype))(idx,:)=[mean(cs1),mean(cs2),mean(es1),mean(es2)];
        idx=idx+1;
    end
    fidx=fidx+1;
end
plot([-10,10],[-10,10],'--k')
xlim([-2,2]);
ylim(xlim());
xlabel('Sample 1 normalized FR (Z-score)')
ylabel('Sample 2 normalized FR (Z-score)')
legend([ch(1),ch(2),eh],{'S1 correct trials','S2 correct trials','Error trials'});


figure('Color','w','Position',[100,100,1200,300]);
subplot(1,3,1);
hold on;
ph=histogram([stats.(sprintf('type%d',types(1)))(:,1);stats.(sprintf('type%d',types(2)))(:,2)],-2:0.2:2,'Normalization','probability','FaceColor','r','FaceAlpha',0.4);
nph=histogram([stats.(sprintf('type%d',types(1)))(:,2);stats.(sprintf('type%d',types(2)))(:,1)],-2:0.2:2,'Normalization','probability','FaceColor','k','FaceAlpha',0.4);
xline(nanmean([stats.(sprintf('type%d',types(1)))(:,1);stats.(sprintf('type%d',types(2)))(:,2)]),'--r','LineWidth',1);
xline(nanmean([stats.(sprintf('type%d',types(1)))(:,2);stats.(sprintf('type%d',types(2)))(:,1)]),'--k','LineWidth',1)
title('Correct trials')
xlabel('Normalized FR (Z-Score)');
ylabel('Probability')



subplot(1,3,2);
hold on;
ph=histogram([stats.(sprintf('type%d',types(1)))(:,3);stats.(sprintf('type%d',types(2)))(:,4)],-2:0.2:2,'Normalization','probability','FaceColor','r','FaceAlpha',0.4);
nph=histogram([stats.(sprintf('type%d',types(1)))(:,4);stats.(sprintf('type%d',types(2)))(:,3)],-2:0.2:2,'Normalization','probability','FaceColor','k','FaceAlpha',0.4);
xline(nanmean([stats.(sprintf('type%d',types(1)))(:,3);stats.(sprintf('type%d',types(2)))(:,4)]),'--r','LineWidth',1)
xline(nanmean([stats.(sprintf('type%d',types(1)))(:,4);stats.(sprintf('type%d',types(2)))(:,3)]),'--k','LineWidth',1)
title('Error trials');
legend([ph,nph],{'Prefered','Non-prefered'});
xlabel('Normalized FR (Z-Score)');
ylabel('Probability')

onecolumn=size(stats.(sprintf('type%d',types(1))),1)+size(stats.(sprintf('type%d',types(2))),1);
[xc,yc,~,aucc]=perfcurve((1:2*onecolumn)>onecolumn,[stats.(sprintf('type%d',types(1)))(:,1);stats.(sprintf('type%d',types(2)))(:,2);stats.(sprintf('type%d',types(1)))(:,2);stats.(sprintf('type%d',types(2)))(:,1)],0);
[xe,ye,~,auce]=perfcurve((1:2*onecolumn)>onecolumn,[stats.(sprintf('type%d',types(1)))(:,3);stats.(sprintf('type%d',types(2)))(:,4);stats.(sprintf('type%d',types(1)))(:,4);stats.(sprintf('type%d',types(2)))(:,3)],0);
subplot(1,3,3);
hold on;
hc=plot(xc,yc,'-r','LineWidth',1);
he=plot(xe,ye,'-k','LineWidth',1);
legend([hc,he],{sprintf('Correct trials AUC=%0.3f',aucc),sprintf('Error trials AUC=%0.3f',auce)})
xlabel('False positive rate (fpr)');
ylabel('True positive rate (tpr)');
return