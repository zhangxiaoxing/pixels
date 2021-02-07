%%overall density selec
%% further classifyed memory unit connections
% conncount=[];
% for bin=1:6
%     load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
%     if bin==1
%         paircount=nnz(diff(pair_reg,1,2));
%     end
%     conncount=[conncount,nnz(diff(reg_chain_S1,1,2)),nnz(diff(reg_chain_S2,1,2))];
%
% end
% mmm=mean(conncount)/paircount;
% stdm=std(conncount/paircount);
% mmm/paircount;

%%overall density nonsel
thres=100;
if false
    

    conncount=nan(114,6,2);
    pairscount=nan(114,6);
    for bin=1:6
        load(sprintf('0813_nonsel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
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
            selS1=conn_chains_S1{bin}(:,1)>ses*100000 & conn_chains_S1{bin}(:,1)<(ses+1)*100000;
            selS2=conn_chains_S2{bin}(:,1)>ses*100000 & conn_chains_S2{bin}(:,1)<(ses+1)*100000;
            conncount(ses,bin,1)=nnz(diff(reg_chains_S1{bin}(selS1,:),1,2));
            conncount(ses,bin,2)=nnz(diff(reg_chains_S2{bin}(selS2,:),1,2));
            pairsel=pair_chains{bin}(:,1)>ses*100000 & pair_chains{bin}(:,1)<(ses+1)*100000;
            pairscount(ses,bin)=nnz(diff(pair_regs{bin}(pairsel,:),1,2));
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
    bzthres=250;
    load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    s1sel=all(bz_spk_count_S1>bzthres,2);
    s2sel=all(bz_spk_count_S2>bzthres,2);
    pair_regs{bin}=pair_reg;
    reg_chains_S1{bin}=bz_conn_reg_S1(s1sel,:);
    reg_chains_S2{bin}=bz_conn_reg_S2(s2sel,:);
    conn_chains_S1{bin}=bz_conn_chain_S1(s1sel,:);
    conn_chains_S2{bin}=bz_conn_chain_S2(s2sel,:);
    pair_chains{bin}=pair_chain;
    pref_chains_S1{bin}= bz_pref_S1(s1sel,:);
    pref_chains_S2{bin}= bz_pref_S2(s2sel,:);
    pref_pairs{bin}=pref_pair;
    clear pair_reg reg_chain_S1 reg_chain_S2 pair_chain ...
        conn_chain_S1 conn_chain_S2 pref_chain_S1 pref_chain_S2 ...
        pref_pair s1sel s2sel
end
end
congru_active_count=nan(114,6,2);
incong_active_count=nan(114,6,2);
congru_inactive_count=nan(114,6,2);
incong_inactive_count=nan(114,6,2);

congru_active_pair=nan(114,6);
incong_active_pair=nan(114,6);
congru_inactive_pair=nan(114,6);
incong_inactive_pair=nan(114,6);

for bin=1:6
    for ses=1:114
        sessel=conn_chains_S1{bin}(:,1)>ses*100000 & conn_chains_S1{bin}(:,1)<(ses+1)*100000;
        congru_active_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2) & pref_chains_S1{bin}(sessel,bin)==pref_chains_S1{bin}(sessel,bin+6) & pref_chains_S1{bin}(sessel,bin)>0);
        incong_active_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2) & pref_chains_S1{bin}(sessel,bin)~=pref_chains_S1{bin}(sessel,bin+6) & pref_chains_S1{bin}(sessel,bin)>0 & pref_chains_S1{bin}(sessel,bin+6)>0);
        congru_inactive_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2) & max(pref_chains_S1{bin}(sessel,1:6),[],2)==max(pref_chains_S1{bin}(sessel,7:12),[],2) & any(pref_chains_S1{bin}(sessel,[bin,bin+6])<1,2));
        incong_inactive_count(ses,bin,1)=nnz(reg_chains_S1{bin}(sessel,1)~=reg_chains_S1{bin}(sessel,2) & max(pref_chains_S1{bin}(sessel,1:6),[],2)~=max(pref_chains_S1{bin}(sessel,7:12),[],2) & any(pref_chains_S1{bin}(sessel,[bin,bin+6])<1,2));
        
        sessel=conn_chains_S2{bin}(:,1)>ses*100000 & conn_chains_S2{bin}(:,1)<(ses+1)*100000;
        congru_active_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2) & pref_chains_S2{bin}(sessel,bin)==pref_chains_S2{bin}(sessel,bin+6) & pref_chains_S2{bin}(sessel,bin)>0);
        incong_active_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2) & pref_chains_S2{bin}(sessel,bin)~=pref_chains_S2{bin}(sessel,bin+6) & pref_chains_S2{bin}(sessel,bin)>0 & pref_chains_S2{bin}(sessel,bin+6)>0);
        congru_inactive_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2) & max(pref_chains_S2{bin}(sessel,1:6),[],2)==max(pref_chains_S2{bin}(sessel,7:12),[],2) & any(pref_chains_S2{bin}(sessel,[bin,bin+6])<1,2));
        incong_inactive_count(ses,bin,2)=nnz(reg_chains_S2{bin}(sessel,1)~=reg_chains_S2{bin}(sessel,2) & max(pref_chains_S2{bin}(sessel,1:6),[],2)~=max(pref_chains_S2{bin}(sessel,7:12),[],2) & any(pref_chains_S2{bin}(sessel,[bin,bin+6])<1,2));
        
        sessel=pair_chains{bin}(:,1)>ses*100000 & pair_chains{bin}(:,1)<(ses+1)*100000;
        congru_active_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & pref_pairs{bin}(sessel,bin)==pref_pairs{bin}(sessel,bin+6) & pref_pairs{bin}(sessel,bin)>0);
        incong_active_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & pref_pairs{bin}(sessel,bin)~=pref_pairs{bin}(sessel,bin+6) & pref_pairs{bin}(sessel,bin)>0 & pref_pairs{bin}(sessel,bin+6)>0);
        congru_inactive_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & max(pref_pairs{bin}(sessel,1:6),[],2)==max(pref_pairs{bin}(sessel,7:12),[],2) & any(pref_pairs{bin}(sessel,[bin,bin+6])<1,2));
        incong_inactive_pair(ses,bin)=nnz(pair_regs{bin}(sessel,1)~=pair_regs{bin}(sessel,2) & max(pref_pairs{bin}(sessel,1:6),[],2)~=max(pref_pairs{bin}(sessel,7:12),[],2) & any(pref_pairs{bin}(sessel,[bin,bin+6])<1,2));
        
    end
end
%%

mca=[];
mia=[];
mcs=[];
mis=[];
for i=1:114
    sescountca=[];
    sescountia=[];
    sescountcs=[];%silent
    sescountis=[];%silent
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
    
end

mmca=mean(mca(:,2));
mmia=mean(mia(:,2));
mmcs=mean(mcs(:,2));
mmis=mean(mis(:,2));

stdca=std(mca(:,2));
stdia=std(mia(:,2));
stdcs=std(mcs(:,2));
stdis=std(mis(:,2));

sum(congru_active_count(congru_active_pair>thres))
sum(congru_active_pair(congru_active_pair>thres))
sum(congru_inactive_count(congru_inactive_pair>thres))
sum(congru_inactive_pair(congru_inactive_pair>thres))
sum(incong_active_count(incong_active_pair>thres))
sum(incong_active_pair(incong_active_pair>thres))
sum(incong_inactive_count(incong_inactive_pair>thres))
sum(incong_inactive_pair(incong_inactive_pair>thres))


% sum(conncount(pairscount>thres))
% sum(pairscount(pairscount>thres))
% 
% mmc=mean(conncount);
% stdc=std(conncount);
% 
% mmi=mean(incongcount);
% stdi=std(incongcount);
% 
% mms=mean(silentcount);
% stds=std(silentcount);

mmn=0;stdn=0;
% close all
fh=figure('Color','w','Position',[100,100,235,270]);
hold on
bar([mmn,mmis,mmia,mmcs,mmca],'FaceColor','w')
errorbar(1:5,[mmn,mmis,mmia,mmcs,mmca],[stdn/sqrt(114),stdis/sqrt(114),stdia/sqrt(114),stdcs/sqrt(114),stdca/sqrt(114)],'k.','CapSize',12)
xlim([0.5,5.5])
ylim([0,0.003]);
set(gca,'XTick',1:5,'XTickLabel',{'Non-memory','Incongruent inactive','Incongruent active','Congruent inactive','Congruent active'},'XTickLabelRotation',30)
ylabel('Connection density')
exportgraphics(gcf,'pairwise-conn-dens_bz.pdf','ContentType','vector');
exportgraphics(fh,'pairwise-conn-dens_bz.png');