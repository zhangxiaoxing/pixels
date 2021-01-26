%% further classified memory unit connections
thres=100;
if true
    conncount=nan(114,6,2);
    pairscount=nan(114,6);
    for bin=1:6
        load(sprintf('0115_nonsel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        pair_regs{bin}=pair_reg;
        reg_chains_S1{bin}=reg_chain_S1;
        reg_chains_S2{bin}=reg_chain_S2;
        conn_chains_S1{bin}=conn_chain_S1;
        conn_chains_S2{bin}=conn_chain_S2;
        pair_chains{bin}=pair_chain;
        clear pair_reg reg_chain_S1 reg_chain_S2 pair_chain
    end
    for bin=1:6
        for ses=1:114
            selS1=conn_chains_S1{bin}(:,1)>=ses*100000 & conn_chains_S1{bin}(:,1)<(ses+1)*100000;
            selS2=conn_chains_S2{bin}(:,1)>=ses*100000 & conn_chains_S2{bin}(:,1)<(ses+1)*100000;
            conncount(ses,bin,1)=nnz(diff(reg_chains_S1{bin}(selS1,:),1,2) & max(reg_chains_S1{bin}(selS1,:),[],2)<116);
            conncount(ses,bin,2)=nnz(diff(reg_chains_S2{bin}(selS2,:),1,2) & max(reg_chains_S2{bin}(selS2,:),[],2)<116);
            pairsel=pair_chains{bin}(:,1)>=ses*100000 & pair_chains{bin}(:,1)<(ses+1)*100000;
            pairscount(ses,bin)=nnz(diff(pair_regs{bin}(pairsel,:),1,2) & max(pair_regs{bin}(pairsel,:),[],2)<116);
        end
    end
    mn=[];
    for i=1:114
        sescountca=[];
        for bin=1:6
            if pairscount(i,bin)>thres
                sescountca=[sescountca;squeeze(conncount(i,bin,:)/pairscount(i,bin))];
            end
        end
        if ~isempty(sescountca)
            mn=[mn;i,mean(sescountca)];
        end
    end
    
    mmn=mean(mn(:,2));
    stdn=std(mn(:,2));
end
%%stats
% [~,~,p]=crosstab((1:244128+1526311)>244128 ,[(1:244128)>46158,(1:1526311)>202185 ])

%% congruent pairs
if true
for bin=1:6
    load(sprintf('0116_memory_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    pair_regs{bin}=pair_reg;
    reg_chains_S1{bin}=reg_chain_S1;
    reg_chains_S2{bin}=reg_chain_S2;
    conn_chains_S1{bin}=conn_chain_S1;
    conn_chains_S2{bin}=conn_chain_S2;
    pair_chains{bin}=pair_chain;
    pref_chains_S1{bin}=pref_chain_S1;
    pref_chains_S2{bin}=pref_chain_S2;
    pref_pairs{bin}=pref_pair;
    clear pair_reg reg_chain_S1 reg_chain_S2 pair_chain conn_chain_S1 conn_chain_S2 pref_chain_S1 pref_chain_S2 pref_pair
end
end
congru_active_count=nan(114,6,2);
incong_active_count=nan(114,6,2);
congru_inactive_count=nan(114,6,2);
incong_inactive_count=nan(114,6,2);
memory_count=nan(114,6,2);

congru_active_pair=nan(114,6);
incong_active_pair=nan(114,6);
congru_inactive_pair=nan(114,6);
incong_inactive_pair=nan(114,6);
memory_pair=nan(114,6);

for bin=1:6
    for ses=1:114
        sessel=conn_chains_S1{bin}(:,1)>ses*100000 & conn_chains_S1{bin}(:,1)<(ses+1)*100000;
        congru_active_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2) & max(reg_chains_S1{bin}(sessel,:),[],2)<116 & pref_chains_S1{bin}(sessel,bin)==pref_chains_S1{bin}(sessel,bin+6) & pref_chains_S1{bin}(sessel,bin)>0);
        incong_active_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2) & max(reg_chains_S1{bin}(sessel,:),[],2)<116 & pref_chains_S1{bin}(sessel,bin)~=pref_chains_S1{bin}(sessel,bin+6) & pref_chains_S1{bin}(sessel,bin)>0 & pref_chains_S1{bin}(sessel,bin+6)>0);
        congru_inactive_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2) & max(reg_chains_S1{bin}(sessel,:),[],2)<116 & max(pref_chains_S1{bin}(sessel,1:6),[],2)==max(pref_chains_S1{bin}(sessel,7:12),[],2) & any(pref_chains_S1{bin}(sessel,[bin,bin+6])<1,2));
        incong_inactive_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2) & max(reg_chains_S1{bin}(sessel,:),[],2)<116 & max(pref_chains_S1{bin}(sessel,1:6),[],2)~=max(pref_chains_S1{bin}(sessel,7:12),[],2) & any(pref_chains_S1{bin}(sessel,[bin,bin+6])<1,2));
        memory_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2)& max(reg_chains_S1{bin}(sessel,:),[],2)<116);
        
        sessel=conn_chains_S2{bin}(:,1)>ses*100000 & conn_chains_S2{bin}(:,1)<(ses+1)*100000;
        congru_active_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2) & max(reg_chains_S2{bin}(sessel,:),[],2)<116 & pref_chains_S2{bin}(sessel,bin)==pref_chains_S2{bin}(sessel,bin+6) & pref_chains_S2{bin}(sessel,bin)>0);
        incong_active_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2) & max(reg_chains_S2{bin}(sessel,:),[],2)<116 & pref_chains_S2{bin}(sessel,bin)~=pref_chains_S2{bin}(sessel,bin+6) & pref_chains_S2{bin}(sessel,bin)>0 & pref_chains_S2{bin}(sessel,bin+6)>0);
        congru_inactive_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2) & max(reg_chains_S2{bin}(sessel,:),[],2)<116 & max(pref_chains_S2{bin}(sessel,1:6),[],2)==max(pref_chains_S2{bin}(sessel,7:12),[],2) & any(pref_chains_S2{bin}(sessel,[bin,bin+6])<1,2));
        incong_inactive_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2) & max(reg_chains_S2{bin}(sessel,:),[],2)<116 & max(pref_chains_S2{bin}(sessel,1:6),[],2)~=max(pref_chains_S2{bin}(sessel,7:12),[],2) & any(pref_chains_S2{bin}(sessel,[bin,bin+6])<1,2));
        memory_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2)& max(reg_chains_S2{bin}(sessel,:),[],2)<116);
        
        sessel=pair_chains{bin}(:,1)>ses*100000 & pair_chains{bin}(:,1)<(ses+1)*100000;
        congru_active_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & max(pair_regs{bin}(sessel,:),[],2)<116 & pref_pairs{bin}(sessel,bin)==pref_pairs{bin}(sessel,bin+6) & pref_pairs{bin}(sessel,bin)>0);
        incong_active_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & max(pair_regs{bin}(sessel,:),[],2)<116 & pref_pairs{bin}(sessel,bin)~=pref_pairs{bin}(sessel,bin+6) & pref_pairs{bin}(sessel,bin)>0 & pref_pairs{bin}(sessel,bin+6)>0);
        congru_inactive_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & max(pair_regs{bin}(sessel,:),[],2)<116 & max(pref_pairs{bin}(sessel,1:6),[],2)==max(pref_pairs{bin}(sessel,7:12),[],2) & any(pref_pairs{bin}(sessel,[bin,bin+6])<1,2));
        incong_inactive_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & max(pair_regs{bin}(sessel,:),[],2)<116 & max(pref_pairs{bin}(sessel,1:6),[],2)~=max(pref_pairs{bin}(sessel,7:12),[],2) & any(pref_pairs{bin}(sessel,[bin,bin+6])<1,2));
        memory_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & max(pair_regs{bin}(sessel,:),[],2)<116);
    end
end
%%

mca=[];
mia=[];
mcs=[];
mis=[];
mem=[];
for i=1:114
    sescountca=[];
    sescountia=[];
    sescountcs=[];%silent
    sescountis=[];%silent
    sescount_mem=[];
    for bin=1:6
        if congru_active_pair(i,bin)>thres
            sescountca=[sescountca;squeeze(congru_active_count(i,bin,:)/congru_active_pair(i,bin))];
        end
        if incong_active_pair(i,bin)>thres
            sescountia=[sescountia;squeeze(incong_active_count(i,bin,:)/incong_active_pair(i,bin))];
        end
        
        if  congru_inactive_pair(i,bin)>thres
            sescountcs=[sescountcs;squeeze(congru_inactive_count(i,bin,:)/congru_inactive_pair(i,bin))];
        end
        if  incong_inactive_pair(i,bin)>thres
            sescountis=[sescountis;squeeze(incong_inactive_count(i,bin,:)/incong_inactive_pair(i,bin))];
        end
        
        
        if memory_pair(i,bin)>thres
            sescount_mem=[sescount_mem;squeeze(memory_count(i,bin,:)/memory_pair(i,bin))];
        end
        
    end
    if ~isempty(sescountca)
        mca=[mca;i,mean(sescountca)];
    end
    
    if ~isempty(sescountia)
        mia=[mia;i,mean(sescountia)];
    end
    
    if ~isempty(sescountcs)
        mcs=[mcs;i,mean(sescountcs)];
    end
    if ~isempty(sescountis)
        mis=[mis;i,mean(sescountis)];
    end
    
    
    if ~isempty(sescount_mem)
        mem=[mem;i,mean(sescount_mem)];
    end
end

mmca=mean(mca(:,2));
mmia=mean(mia(:,2));
mmcs=mean(mcs(:,2));
mmis=mean(mis(:,2));
mmem=mean(mem(:,2));

stdca=std(mca(:,2));
stdia=std(mia(:,2));
stdcs=std(mcs(:,2));
stdis=std(mis(:,2));

sum(congru_active_count(cat(3,congru_active_pair,congru_active_pair)>thres))
sum(congru_active_pair(congru_active_pair>thres)).*2
sum(congru_inactive_count(cat(3,congru_inactive_pair,congru_inactive_pair)>thres))
sum(congru_inactive_pair(congru_inactive_pair>thres)).*2
sum(incong_active_count(cat(3,incong_active_pair,incong_active_pair)>thres))
sum(incong_active_pair(incong_active_pair>thres)).*2
sum(incong_inactive_count(cat(3,incong_inactive_pair,incong_inactive_pair)>thres))
sum(incong_inactive_pair(incong_inactive_pair>thres)).*2


sum(conncount(cat(3,pairscount,pairscount)>thres))
sum(pairscount(pairscount>thres)).*2

sum(conncount(cat(3,pairscount,pairscount)>thres))
sum(pairscount(pairscount>thres)).*2

%%ANOVA test for 1K
keyboard()
cnt={conncount,incong_inactive_count,incong_active_count,congru_inactive_count,congru_active_count};
pairs={pairscount,incong_inactive_pair,incong_active_pair,congru_inactive_pair,congru_active_pair};

y=[];
g=[];
for i=1:5
    for ss=1:2
        cc=cnt{i}(:,:,ss);
        pp=pairs{i};
        sel=pp>thres;
        y=[y;cc(sel)./pp(sel)];
        g=[g;i.*ones(nnz(sel),1)];
    end
end
anova1(y,g)

%%chisq test for text
keyboard()
sum(memory_count(cat(3,memory_pair,memory_pair)>thres))
sum(congru_active_pair(congru_active_pair>thres)).*2


% 
% mmc=mean(conncount);
% stdc=std(conncount);
% 
% mmi=mean(incongcount);
% stdi=std(incongcount);
% 
% mms=mean(silentcount);
% stds=std(silentcount);


% close all
fh=figure('Color','w','Position',[100,100,235,270]);
hold on
bar([mmn,mmis,mmia,mmcs,mmca].*100,'FaceColor','w')
errorbar(1:5,[mmn,mmis,mmia,mmcs,mmca].*100,[stdn/sqrt(114),stdis/sqrt(114),stdia/sqrt(114),stdcs/sqrt(114),stdca/sqrt(114)].*100,'k.','CapSize',12)
xlim([0.5,5.5])
ylim([0,25]);
set(gca,'XTick',1:5,'XTickLabel',{'Non-memory','Incongruent inactive','Incongruent active','Congruent inactive','Congruent active'},'XTickLabelRotation',30)
ylabel('Functional coupling density (%)')
exportgraphics(gcf,'pairwise-conn-dens.pdf','ContentType','vector');
exportgraphics(fh,'pairwise-conn-dens.png');



onebin=pairscount(:,1);
totalpairs=sum(onebin(onebin>thres));

onebin=memory_pair(:,1);
memory_pairs=sum(onebin(onebin>thres));
