to_plot=false;
plot_sel=false;

close all
% load('join_reg_set.mat')
load('reg_keep.mat')
load('zcmap.mat')

% 
% allkeep=~ismember(join_reg_set,{'Unlabeled','root'});
% all_reg=join_reg_set(allkeep);
all_reg=reg_set;
major=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), all_reg);

conn_list=cell(6,1);
pair_list=cell(6,1);
c_mat_sum=zeros(length(all_reg));
p_mat_sum=zeros(length(all_reg));
zero_count=zeros(size(all_reg));
for bin=1:6
    load(sprintf('0626_selec_conn_mat_duo_6s_%d_%d.mat',bin,bin+1))
    load(sprintf('0626_selec_pair_mat_duo_6s_%d_%d.mat',bin,bin+1))
    conn_list{bin}=conn_mat_S1;
    pair_list{bin}=pair_mat;
    for i=1:length(all_reg)
        zero_count(i)=zero_count(i)+nnz(pair_mat(i,:)<20)+nnz(pair_mat(:,i)<20);
    end
    c_mat_sum=c_mat_sum+conn_mat_S1;
    p_mat_sum=p_mat_sum+pair_mat;
    
end

% ckeep=(zero_count<1500 & major);
% reg_keep=all_reg(ckeep);

p_mat_sum(p_mat_sum<50)=0;
nz_p_count=arrayfun(@(x) sum(nnz([p_mat_sum(:,x);p_mat_sum(x,:)'])),1:length(p_mat_sum));
[~,I] = sort(nz_p_count(:),'descend');
ckeep=I(1:30);
reg_keep=all_reg(ckeep);

ratio_mat=c_mat_sum(ckeep,ckeep)./p_mat_sum(ckeep,ckeep);
% rsum=arrayfun(@(x) nanmean([ratio_mat(:,x);ratio_mat(x,[1:x-1,x+1:end])']),1:length(ratio_mat));
% [~,leafOrder]=sort(rsum,'descend');

D=pdist(ratio_mat',@nan_spearman);
Z=linkage(D,'average');
leafOrder = optimalleaforder(Z,D,'Criteria','group');
% figure('Color','w')
% [H,T,outperm]=dendrogram(Z,4,'Orientation','left','ColorThreshold','default','Reorder',leafOrder,'Labels',keep_reg);

% fh=figure('Color','w','Position',[100,100,420,380]);
% ratio_matfig=ratio_mat(leafOrder,leafOrder);
% 
% h=imagesc(ratio_matfig,[0,0.3]);
% set(h,'alphadata',~isnan(ratio_matfig)); 
% colormap(zcmap./255);
% colorbar
% set(gca,'FontSize',6.5,'xaxisLocation','top','XTick',1:nnz(ckeep),'YTick',1:nnz(ckeep),'XTickLabel',reg_keep(leafOrder),'YTickLabel',reg_keep(leafOrder),'XTickLabelRotation',90,'Color',[58 56 56]./255);
% keyboard
% exportgraphics(fh,sprintf('ratio_mat_6s_%d_%d.pdf',bin,bin+1),'Resolution',300)


ratio_list=cell(1,6);
ratio_hist=[];
if plot_sel
    for bin=1:6
        conn_mat=conn_list{bin};
        pair_mat=pair_list{bin};
        pair_mat(pair_mat<20)=0;
        conn_mat(pair_mat<20)=0;
%         keyboard
        ratio_mat=conn_mat./pair_mat;
        ratio_hist=[ratio_hist,ratio_mat(:)];
        
        if to_plot
            fh=figure('Color','w','Position',[100,100,420,380]);
            ratio_mat=ratio_mat(ckeep,ckeep);
            h=imagesc(ratio_mat(leafOrder,leafOrder),[0,0.25]);
            set(h,'alphadata',~isnan(ratio_mat)); 
            colormap(zcmap./255);
            colorbar
            set(gca,'FontSize',6.5,'xaxisLocation','top','XTick',1:nnz(ckeep),'YTick',1:nnz(ckeep),'XTickLabel',reg_keep(leafOrder),'YTickLabel',reg_keep(leafOrder),'XTickLabelRotation',90,'Color',[58 56 56]./255);
            exportgraphics(fh,sprintf('ratio_mat_6s_%d_%d.pdf',bin,bin+1),'Resolution',300)
%             print(fh,sprintf('ratio_mat_6s_%d_%d.png',bin,bin+1),'-dpng','-r300')
        end
        
        
        pair_mat(pair_mat<500)=0;
        conn_mat(pair_mat<500)=0;
        ratio_list{bin}=conn_mat./pair_mat;
    end
end
close all
fh=figure('Color','w','Position',[100,100,235,130])
histogram(ratio_hist,0:0.02:0.5,'Normalization','probability','FaceColor','k');
ylim([0,0.02])
xlabel('connection density')
ylabel('normalized probability')
set(gca,'YTick',0:0.01:0.02)
exportgraphics(fh,'conn_dense_hist.pdf','ContentType','vector')


nonsel_ratio_list=cell(1,6);
plot_non_sel=false;
if plot_non_sel
    for bin=1
        nsconnfstr=load(sprintf('nonsel_conn_mat_duo_6s_%d_%d.mat',bin,bin+1));
        nspairfstr=load(sprintf('nonsel_pair_mat_duo_6s_%d_%d.mat',bin,bin+1));
        conn_mat=nsconnfstr.conn_mat;
        pair_mat=nspairfstr.pair_mat;
        ratio_mat=(conn_mat(ckeep,ckeep)./pair_mat(ckeep,ckeep)+minv/2);
        nonsel_ratio_list{bin}=ratio_mat(leafOrder,leafOrder);
        if to_plot
            fh=figure('Color','w','Position',[100,100,420,380]);
            h=imagesc(ratio_mat(leafOrder,leafOrder),[-3,0]);
            set(h,'alphadata',~isnan(ratio_mat)); 
            colormap(zcmap./255);
            colorbar
            set(gca,'FontSize',6.5,'xaxisLocation','top','XTick',1:nnz(ckeep),'YTick',1:nnz(ckeep),'XTickLabel',keep_reg(leafOrder),'YTickLabel',keep_reg(leafOrder),'XTickLabelRotation',90,'Color',[58 56 56]./255);
            exportgraphics(fh,sprintf('nonsel_ratio_mat_6s_%d_%d.pdf',bin,bin+1),'Resolution',300)
%             print(fh,sprintf('nonsel_ratio_mat_6s_%d_%d.png',bin,bin+1),'-dpng','-r300')
        end
    end
end

sel_corr_bin=zeros(5,1);
r_list=nan(6,6);
for bin=1:5
    for binnext=(bin+1):6
        binmat=ratio_list{bin};
        nextmat=ratio_list{binnext};
        corrmat=zeros(0,2);
        for i=1:length(binmat)
            for j=1:length(binmat)
                if ~isnan(binmat(i,j)) && ~isnan(nextmat(i,j))
                    corrmat(end+1,:)=[binmat(i,j),nextmat(i,j)];
                end
            end
        end
        [r,p]=corrcoef(corrmat(:,1),corrmat(:,2));
        r_list(bin,binnext)=r(1,2);
    end
end
fh=figure('Color','w','Position',[100,100,480,100]);
h=imagesc(r_list,[0,1]);
colormap('jet');
colorbar();
set(gca(),'YDir','normal','Color','w')
set(h,'alphadata',~isnan(r_list)); 
xlabel('bin #')
ylabel('bin #')
exportgraphics(fh,sprintf('bin_bin_corr_6s.pdf',bin,bin+1),'Resolution',300)

return
%% plot sel nonsel corr coef
for bin=1
    sel_mat=ratio_list{bin};
    nonsel_mat=nonsel_ratio_list{bin};
    corrmat=zeros(0,2);
    for i=1:length(sel_mat)
        for j=1:length(sel_mat)
            if ~isnan(sel_mat(i,j)) && ~isnan(nonsel_mat(i,j))
                corrmat(end+1,:)=[sel_mat(i,j),nonsel_mat(i,j)];
            end
        end
    end
    
    fh=figure('Color','w','Position',[100,100,250,250]);
    hold on;
    ph=scatter(corrmat(:,1),corrmat(:,2),15,'MarkerEdgeColor','none','MarkerFaceAlpha',0.25,'MarkerFaceColor','k');
    plot([-5,0],[-5,0],'--r');
    ax=gca();
    xlim([-3.5,0])
    ylim([-3.5,0])
    xlabel('selective neuron norm. connection density');
    ylabel('nonselective neuron norm. connection density');
    ax.XTick=-3:0;
    ax.YTick=-3:0;
    [r,p]=corrcoef(corrmat(:,1),corrmat(:,2));
    exportgraphics(fh,sprintf('sel_nonsel_corr_6s_%d_%d.pdf',bin,bin+1),'Resolution',300)
end




return

for bin=1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Generate data export for gephi network figure%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conn_mat=conn_list{bin};
    pair_mat=pair_list{bin};
    ratio_mat=conn_mat(ckeep,ckeep)./pair_mat(ckeep,ckeep);
    
    csvcell=cell(1,4);
    csvcell(1,:)={'Id','Label','Interval','Partition'};
    for i=1:length(ratio_mat)
        csvcell(end+1,:)={keep_reg{i},keep_reg{i},[],T(i)};
    end
    
    writecell(csvcell,sprintf('node_exp_%d.csv',bin));
    
    csvcell=cell(1,3);
    csvcell(1,:)={'Source','Target','Weight'};
    for i=1:length(ratio_mat)
        for j=1:length(ratio_mat)
            if ~isnan(ratio_mat(i,j)) && ~(i==j)
                csvcell(end+1,:)={keep_reg{j},keep_reg{i},ratio_mat(i,j)};
            end
        end
    end
    
    writecell(csvcell,sprintf('ratio_mat_conn_%d.csv',bin));
    %%%%%%%%%%%%%% End of GEPHI EXPORT %%%%%%%%%%%%%%%%%%%%%%%
%     keyboard
    
%     G = digraph(ratio_mat,keep_reg,'omitselfloops');
%     fh=figure('Color','w');
%     plot(G, 'LineWidth', G.Edges.Weight*5,'EdgeAlpha',0.5,'ArrowSize',6,'Layout','force')
%     keyboard
%     set(fh,'visible','off')
%     set(fh,'PaperSize',[15,10])
%     print(fh,sprintf('connection_force_%d_%d.pdf',bin_range(1),bin_range(2)),'-dpdf','-r300')
%     set(fh, 'PaperPosition', [0 0 12 12])
%     print(fh,sprintf('connection_force_%d_%d.png',bin_range(1),bin_range(2)),'-dpng','-r300')

end


function  out=nan_spearman(XI,XJ)
if min(size(XJ))==1
    sel= ~(isnan(XI) | isnan(XJ)| isinf(XI) | isinf(XJ));
    if nnz(sel)<3
        out=2;
    else
        NI=XI(sel);
        NJ=XJ(sel);
        nand=corr(NI(:),NJ(:),'type','Spearman');
        out=1-nand;
        if isnan(out)
            out=2;
        end
    end

else
    nand=zeros(size(XJ,1),1);
    for i=1:size(XJ,1)
%         disp(i)
        sel= ~(isnan(XI) | isnan(XJ(i,:))| isinf(XI) | isinf(XJ(i,:)));
        if nnz(sel)<3
            nand(i)=-1;
        else
            NI=XI(sel);
            NJ=XJ(i,sel);
            roh=corr(NI(:),NJ(:),'type','Spearman');
            nand(i)=roh;
            if isnan(roh)
                nand(i)=-1;
            end
        end
    end
    out=1-nand;
end
end


