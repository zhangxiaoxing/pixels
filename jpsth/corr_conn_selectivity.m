fstr=cell(1,6);
conn_all=[];
congru_all=[];
pref_all=[];
for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    t=fstr{bin}.conn_chain_S1(:,[1,2,2,2]);
    t(:,3)=1;
    t(:,4)=bin;
    conn_all=[conn_all;t];
    sel=fstr{bin}.pref_chain_S1(:,bin)==fstr{bin}.pref_chain_S1(:,bin+6) & fstr{bin}.pref_chain_S1(:,bin)>0;
    congru_all=[congru_all;t(sel,:)];
    pref_one=[fstr{bin}.pair_chain(:,1),fstr{bin}.pref_pair(:,1:6);...
        fstr{bin}.pair_chain(:,2),fstr{bin}.pref_pair(:,7:12)];
    pref_all=[pref_all;pref_one];
    
    t=fstr{bin}.conn_chain_S2(:,[1,2,2,2]);
    t(:,3)=2;
    t(:,4)=bin;
    conn_all=[conn_all;t];
    sel=fstr{bin}.pref_chain_S2(:,bin)==fstr{bin}.pref_chain_S2(:,bin+6) & fstr{bin}.pref_chain_S2(:,bin)>0;
    congru_all=[congru_all;t(sel,:)];
%     pref_one=[fstr{bin}.conn_chain_S2(:,1),fstr{bin}.pref_chain_S2(:,1:6);...
%         fstr{bin}.conn_chain_S2(:,2),fstr{bin}.pref_chain_S2(:,7:12)];
%     pref_all=[pref_all;pref_one];
end
pref_all=unique(pref_all,'rows');
%% all connection
corr_list=[];
for idx=1:length(pref_all)
    if rem(idx,100)==0
        disp(idx)
    end
    one=[pref_all(idx,1),nnz(pref_all(idx,2:end)),nnz(conn_all(:,2)==pref_all(idx,1))];
    corr_list=[corr_list;one];
end

fh=figure('Color','w');
hold on
lh=scatter(corr_list(:,2)+rand(length(corr_list),1)*0.5-0.25,corr_list(:,3),3,'.');
for i=1:6
    sel=corr_list(:,2)==i;
    plot(i,mean(corr_list(sel,3)),'ro')
end
xlim([0,7])
xlabel('Selective bins')
ylabel('Total pre-unit * bins')

%% congru connection
congru_list=[];
for idx=1:length(pref_all)
    if rem(idx,100)==0
        disp(idx)
    end
    one=[pref_all(idx,1),nnz(pref_all(idx,2:end)),nnz(congru_all(:,2)==pref_all(idx,1))];
    congru_list=[congru_list;one];
end

fh=figure('Color','w');
hold on
scatter(congru_list(:,2)+rand(length(congru_list),1)*0.5-0.25,congru_list(:,3),3,'.');
for i=1:6
    sel=congru_list(:,2)==i;
    plot(i,mean(congru_list(sel,3)),'ro')
end
xlim([0,7])
xlabel('Selective bins')
ylabel('Congruent pre-unit * bins')


%% congru connection
congru_list=[];
for idx=1:length(pref_all)
    if rem(idx,100)==0
        disp(idx)
    end
    one=[pref_all(idx,1),nnz(pref_all(idx,2:end)),nnz(unique(congru_all(congru_all(:,2)==pref_all(idx,1),1)))];
    congru_list=[congru_list;one];
end

fh=figure('Color','w');
hold on
for i=1:6
    sel=congru_list(:,2)==i;
    mm=mean(congru_list(sel,3));
    ci=bootci(1000,@(x) mean(x),congru_list(sel,3));
    bar(i,mm,'FaceColor','w');
    errorbar(i,mm,ci(1)-mm,ci(2)-mm,'k.')
end
xlim([0.5,6.5])
xlabel('Memory neurons with selective bins')
ylabel('Number of unique congruent pre-units')
fh.Position(3:4)=[235,235];
set(gca,'YTick',0:10:20,'XTick',1:6);
exportgraphics(fh,'corr_sel_congru_receive.pdf')




%% congru post
congru_list=[];
for idx=1:length(pref_all)
    if rem(idx,100)==0
        disp(idx)
    end
    one=[pref_all(idx,1),nnz(pref_all(idx,2:end)),nnz(unique(congru_all(congru_all(:,1)==pref_all(idx,1),2)))];
    congru_list=[congru_list;one];
end

fh=figure('Color','w');
hold on
for i=1:6
    sel=congru_list(:,2)==i;
    mm=mean(congru_list(sel,3));
    ci=bootci(1000,@(x) mean(x),congru_list(sel,3));
    bar(i,mm,'FaceColor','w');
    errorbar(i,mm,ci(1)-mm,ci(2)-mm,'k.')
end
xlim([0.5,6.5])
xlabel('Memory neurons with selective bins')
ylabel('Number of unique congruent post-units')
fh.Position(3:4)=[235,235];
set(gca,'YTick',0:10:20,'XTick',1:6);
exportgraphics(fh,'corr_sel_congru_output.pdf')





%% rings

load rings.mat
allrings=cell2mat(rings(:));
ring_list=[];
for idx=1:length(pref_all)
    if rem(idx,100)==0
        disp(idx)
    end
    one=[pref_all(idx,1),nnz(pref_all(idx,2:end)),nnz(allrings(:)==pref_all(idx,1))];
    ring_list=[ring_list;one];
end

fh=figure('Color','w');
hold on
% scatter(ring_list(:,2)+rand(length(ring_list),1)*0.5-0.25,ring_list(:,3),3,'.');
for i=1:6
    sel=ring_list(:,2)==i;
    mm=mean(ring_list(sel,3));
    ci=bootci(1000,@(x) mean(x),ring_list(sel,3));
    bar(i,mm,'FaceColor','w');
    errorbar(i,mm,ci(1)-mm,ci(2)-mm,'k.')
end
ylim([0,8])
xlim([0.5,6.5])
xlabel('Memory neurons with selective bins')
ylabel('Apperance in unique triplet rings')
fh.Position(3:4)=[235,235];
set(gca,'YTick',[0 5]);
exportgraphics(fh,'corr_sel_rings.pdf')



%% rings4
load rings.mat
allrings=cell2mat(rings4(:));
allrings=allrings(:,1:4);
% nnz(arrayfun(@(x) numel(unique(allrings(x,:))),1:length(allrings))~=4)
ring_list=[];
for idx=1:length(pref_all)
    if rem(idx,100)==0
        disp(idx)
    end
    one=[pref_all(idx,1),nnz(pref_all(idx,2:end)),nnz(allrings(:)==pref_all(idx,1))];
    ring_list=[ring_list;one];
end

fh=figure('Color','w');
hold on
% scatter(ring_list(:,2)+rand(length(ring_list),1)*0.5-0.25,ring_list(:,3),3,'.');
for i=1:6
    sel=ring_list(:,2)==i;
    mm=mean(ring_list(sel,3));
    ci=bootci(1000,@(x) mean(x),ring_list(sel,3));
    bar(i,mm,'FaceColor','w');
    errorbar(i,mm,ci(1)-mm,ci(2)-mm,'k.')
end
ylim([0,60])
xlim([0.5,6.5])
xlabel('Memory neurons with selective bins')
ylabel('Apperance in unique quadruplet rings')
fh.Position(3:4)=[235,235];
set(gca,'YTick',[0,50]);
exportgraphics(fh,'corr_sel_rings4.pdf')

