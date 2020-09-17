use_nonsel=false;
binRange=6;


load('zcmap.mat')

greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
greg=reg_set(greymatter);

%% memory, nonmemory-diff
[connmapM,pairmapM]=plot_connmap(false,6);
[connmapN,pairmapN]=plot_connmap(true,4);

%chisq
csvcell=cell(0,8);
rpt_measure=0;
for i=1:length(connmapM)
    for j=1:length(connmapN)
        if i==j || pairmapM(i,j)<100 || pairmapN(i,j)<100
            continue
        end
        rpt_measure=rpt_measure+1;
        [~,~,p]=crosstab((1:(pairmapM(i,j)+pairmapN(i,j)))>pairmapM(i,j),[(1:pairmapM(i,j))>connmapM(i,j),(1:pairmapN(i,j))>connmapN(i,j)]);
        if p*rpt_measure<0.05
            diffIdx=(connmapM(i,j)/pairmapM(i,j)-connmapN(i,j)/pairmapN(i,j))/(connmapM(i,j)/pairmapM(i,j)+connmapN(i,j)/pairmapN(i,j));
            csvcell(end+1,:)={greg{j},greg{i},connmapM(i,j)/pairmapM(i,j),diffIdx,connmapM(i,j),pairmapM(i,j),connmapN(i,j),pairmapN(i,j)};
        end
    end
end
[~,idx]=sort(abs(cell2mat(csvcell(:,4))),'descend');
csvcell=csvcell(idx,:);
writecell(csvcell(:,1:4),'selec_nonsel_diff_idx.csv');
%TODO export



%% reciprocal corr
if false
    [connmap,pairmap]=plot_connmap(use_nonsel,binRange,greg);
    plot_recip(connmap,pairmap)
    keyboard
end
to_plot=false;
plot_hist=true;

function [c_mat_sum,p_mat_sum]= plot_connmap(use_nonsel,binRange,greg)

close all


conn_list=cell(6,1);
pair_list=cell(6,1);
c_mat_sum=zeros(length(greg));
p_mat_sum=zeros(length(greg));
zero_count=zeros(size(greg));
for bin=1:binRange %%FIXME: needs to extends window to full 6s
    if use_nonsel
        load(sprintf('0729_nonsel_conn_mat_duo_6s_%d_%d.mat',bin,bin+1))
        load(sprintf('0729_nonsel_pair_mat_duo_6s_%d_%d.mat',bin,bin+1))
    else
        load(sprintf('0626_selec_conn_mat_duo_6s_%d_%d.mat',bin,bin+1))
        load(sprintf('0810_selec_pair_mat_duo_6s_%d_%d.mat',bin,bin+1))
    end
    conn_mat_S1=conn_mat_S1(greymatter,greymatter);
    conn_list{bin}=conn_mat_S1;
    pair_mat=pair_mat(greymatter,greymatter);
    pair_list{bin}=pair_mat;
    for i=1:length(greg)
        zero_count(i)=zero_count(i)+nnz(pair_mat(i,:)<20)+nnz(pair_mat(:,i)<20);
    end
    c_mat_sum=c_mat_sum+conn_mat_S1;
    p_mat_sum=p_mat_sum+pair_mat;
    
end

%% use sum data for figure plot
p_mat_sum(p_mat_sum<100)=0;
c_mat_sum(p_mat_sum<100)=0;
%%
if true
    export_data_csv(c_mat_sum,p_mat_sum,greg)
    return
end

if use_nonsel
    leafOrder=evalin('base','leafOrder');
    ckeep=evalin('base','ckeep');
    reg_keep=evalin('base','reg_keep');
    ratio_mat=c_mat_sum(ckeep,ckeep)./p_mat_sum(ckeep,ckeep);
else
    nz_p_count=arrayfun(@(x) sum(nnz([p_mat_sum(:,x);p_mat_sum(x,:)'])),1:length(p_mat_sum));
    [~,I] = sort(nz_p_count(:),'descend');
    ckeep=I(1:32);
    reg_keep=greg(ckeep);
    ratio_mat=c_mat_sum(ckeep,ckeep)./p_mat_sum(ckeep,ckeep);
    
    rsum=arrayfun(@(x) nanmean([ratio_mat(:,x);ratio_mat(x,[1:x-1,x+1:end])']),1:length(ratio_mat));
    [~,leafOrder]=sort(rsum,'descend');
    assignin('base','leafOrder',leafOrder);
    assignin('base','ckeep',ckeep);
    assignin('base','reg_keep',reg_keep);
end


% D=pdist(ratio_mat',@nan_spearman);
% Z=linkage(D,'average');
% leafOrder = optimalleaforder(Z,D,'Criteria','group');
% figure('Color','w')
% [H,T,outperm]=dendrogram(Z,4,'Orientation','left','ColorThreshold','default','Reorder',leafOrder,'Labels',keep_reg);

fh=figure('Color','w','Position',[100,100,420,380]);
ratio_matfig=ratio_mat(leafOrder,leafOrder);

h=imagesc(ratio_matfig,[0,0.165]);
set(h,'alphadata',~isnan(ratio_matfig)); 
colormap(zcmap./255);
cbh=colorbar;
cbh.Ticks=0:0.05:0.15;
set(gca,'FontSize',6.5,'xaxisLocation','top','XTick',1:nnz(ckeep),'YTick',1:nnz(ckeep),'XTickLabel',reg_keep(leafOrder),'YTickLabel',reg_keep(leafOrder),'XTickLabelRotation',90,'Color',[58 56 56]./255);
% keyboard
if use_nonsel
    exportgraphics(fh,sprintf('nonsel_avg_ratio_mat_6s_%d_%d.pdf',1,7),'Resolution',300)
else
    exportgraphics(fh,sprintf('selec_sum_avg_ratio_mat_6s_%d_%d.pdf',1,7),'Resolution',300)
end
end

function prob_hist()
ratio_list=cell(1,6);
ratio_hist=[];
if plot_hist
    for bin=1:numel(conn_list)
        conn_mat=conn_list{bin};
        pair_mat=pair_list{bin};
        pair_mat(pair_mat<100)=0;
        conn_mat(pair_mat<100)=0;
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

%% figure 2D prob hist
close all
fh=figure('Color','w','Position',[100,100,235,130]);
ratio_hist=ratio_hist(~isnan(ratio_hist));
histogram(ratio_hist,0:0.01:0.3,'Normalization','probability','FaceColor','k');

xlim([0,0.3])
xlabel('connection density')
ylabel('normalized probability')
ylim([0,0.15])
set(gca,'YTick',0:0.1:0.2)

if use_nonsel
    exportgraphics(fh,'nonsel_conn_dense_hist.pdf','ContentType','vector')
else
    exportgraphics(fh,'selec_conn_dense_hist.pdf','ContentType','vector')
end
return
end
% nonsel_ratio_list=cell(1,6);
% plot_non_sel=false;
% if plot_non_sel
%     for bin=1
%         nsconnfstr=load(sprintf('nonsel_conn_mat_duo_6s_%d_%d.mat',bin,bin+1));
%         nspairfstr=load(sprintf('nonsel_pair_mat_duo_6s_%d_%d.mat',bin,bin+1));
%         conn_mat=nsconnfstr.conn_mat;
%         pair_mat=nspairfstr.pair_mat;
%         ratio_mat=(conn_mat(ckeep,ckeep)./pair_mat(ckeep,ckeep)+minv/2);
%         nonsel_ratio_list{bin}=ratio_mat(leafOrder,leafOrder);
%         if to_plot
%             fh=figure('Color','w','Position',[100,100,420,380]);
%             h=imagesc(ratio_mat(leafOrder,leafOrder),[-3,0]);
%             set(h,'alphadata',~isnan(ratio_mat)); 
%             colormap(zcmap./255);
%             colorbar
%             set(gca,'FontSize',6.5,'xaxisLocation','top','XTick',1:nnz(ckeep),'YTick',1:nnz(ckeep),'XTickLabel',keep_reg(leafOrder),'YTickLabel',keep_reg(leafOrder),'XTickLabelRotation',90,'Color',[58 56 56]./255);
%             exportgraphics(fh,sprintf('nonsel_ratio_mat_6s_%d_%d.pdf',bin,bin+1),'Resolution',300)
% %             print(fh,sprintf('nonsel_ratio_mat_6s_%d_%d.png',bin,bin+1),'-dpng','-r300')
%         end
%     end
% end
function other_stats()
sel_corr_bin=zeros(5,1);
r_list=nan(6,6);
for bin=1:6
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
end

function export_data_csv(c_mat_sum,p_mat_sum,greg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Generate data export for gephi network figure%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assume p_mat_sum(p_mat_sum<100)=0; c_mat_sum(p_mat_sum<100)=0; ckeep; in memory

ratio_mat=c_mat_sum./p_mat_sum;
csvcell=cell(0,3);

%     csvcell(1,:)={'Source','Target','Weight'};
for i=1:length(ratio_mat)
    for j=1:length(ratio_mat)
        if ~isnan(ratio_mat(i,j)) && ~(i==j)
            csvcell(end+1,:)={greg{j},greg{i},ratio_mat(i,j)};
        end
    end
end
[~,idx]=sort(cell2mat(csvcell(:,3)),'descend');
csvcell=csvcell(idx,:);
writecell(csvcell,'selec_sum_ratio.csv');
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

function plot_recip(connmap,pairmap)

plt_list=[];
for i=1:length(connmap)
    for j=1:length(pairmap)
        if i~=j && pairmap(i,j)>=100
            plt_list=[plt_list;connmap(i,j)/pairmap(i,j),connmap(j,i)/pairmap(i,j)];
        end
    end
end
figure('Color','w','Position',[100,100,230,200])
hold on;
scatter(plt_list(:,1),plt_list(:,2),3,'MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
plot([0,0.25],[0,0.25],'r:');
xlim([0,0.25]);
ylim([0,0.25]);
xlabel('in degree')
ylabel('out degree')
set(gca,'XTick',0:0.1:0.2,'YTick',0:0.1:0.2);
exportgraphics(gcf,'in_out_deg_corr.pdf');

[r,p]=corr(plt_list(:,1),plt_list(:,2))
end
