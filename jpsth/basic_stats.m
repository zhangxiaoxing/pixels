%%overall density selec

conncount=[];
for bin=1:6
    load(sprintf('0810_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    if bin==1
        paircount=nnz(diff(pair_reg,1,2));
    end
    conncount=[conncount,nnz(diff(reg_chain_S1,1,2)),nnz(diff(reg_chain_S2,1,2))];
    
end
mmm=mean(conncount)/paircount;
stdm=std(conncount/paircount);
% mmm/paircount;

%%overall density nonsel

conncount=[];
for bin=1:4 %%FIXME: expand to all 6 bins
    load(sprintf('0813_nonsel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    if bin==1
        paircount=nnz(diff(pair_reg,1,2));
    end
    conncount=[conncount,nnz(diff(reg_chain_S1,1,2)),nnz(diff(reg_chain_S2,1,2))];
%     conncount=[conncount,length(reg_chain_S1),length(reg_chain_S2)];
end
mmn=mean(conncount)/paircount;
stdn=std(conncount/paircount);
% mmn/paircount;

%%stats
[~,~,p]=crosstab((1:244128+1526311)>244128 ,[(1:244128)>46158,(1:1526311)>202185 ])

%% congruent pairs
conncount=[];
for bin=1:6
    load(sprintf('0814_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
%     if bin==1
%         paircount=nnz(diff(pair_reg,1,2));
%     end
    congruS1=nnz(reg_chain_S1(:,1)~=reg_chain_S1(:,2) & pref_chain_S1(:,bin)==pref_chain_S1(:,bin+6) & pref_chain_S1(:,bin)>0);
    congruS2=nnz(reg_chain_S2(:,1)~=reg_chain_S2(:,2) & pref_chain_S2(:,bin)==pref_chain_S2(:,bin+6) & pref_chain_S2(:,bin)>0);
    congruPair=nnz(pair_reg(:,1)~=pair_reg(:,2) & pref_pair(:,bin)==pref_pair(:,bin+6) & pref_pair(:,bin)>0);
    conncount=[conncount,congruS1/congruPair,congruS2/congruPair];
    
end
mmc=mean(conncount);
stdc=std(conncount);

close all
figure('Color','w','Position',[100,100,230,130])
hold on
bar([mmn,mmm,mmc],'FaceColor','k')
errorbar(1:3,[mmn,mmm,mmc],[stdn/sqrt(12),stdm/sqrt(12),stdc/sqrt(12)],'k.','CapSize',16,'Color',[0.5,0.5,0.5])
xlim([0.5,3.5])
ylim([0,0.25]);
set(gca,'XTick',1:3,'XTickLabel',{'non-memory','memory units','congruent units'},'XTickLabelRotation',15)
ylabel('connection density')
exportgraphics(gcf,'pairwise-conn-dens.pdf','ContentType','vector');