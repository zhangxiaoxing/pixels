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


function gen_conn_map()
%% post stats statistics
load('reg_keep.mat');
reg_set=reg_set(1:115);

for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
    gen_pair_mat=true;
    if gen_pair_mat
        pair_mat=zeros(length(reg_set),length(reg_set),6);
        for bin=1:6
            for i=1:length(pair_mat)
                disp(i);
                for j=1:length(pair_mat)
                    pair_mat(i,j,bin)=nnz(fstr{bin}.pair_reg(:,1)==i & fstr{bin}.pair_reg(:,2)==j);
                end
            end
        end
%         save('0831_pair_mat_6s.mat','pair_mat');
    end
    gen_conn_mat=true;
    
    if gen_conn_mat
        conn_mat_S1=zeros(length(reg_set),length(reg_set),6);
        conn_mat_S2=zeros(length(reg_set),length(reg_set),6);
        for bin=1:6
            for i=1:length(conn_mat_S1)
                disp(i);
                for j=1:length(conn_mat_S1)
                    conn_mat_S1(i,j,bin)=nnz(fstr{bin}.reg_chain_S1(:,1)==i & fstr{bin}.reg_chain_S1(:,2)==j ...
                        & fstr{bin}.pref_chain_S1(:,bin)==fstr{bin}.pref_chain_S1(:,bin+6) & fstr{bin}.pref_chain_S1(:,bin)>0);
                    conn_mat_S2(i,j,bin)=nnz(fstr{bin}.reg_chain_S2(:,1)==i & fstr{bin}.reg_chain_S2(:,2)==j ...
                        & fstr{bin}.pref_chain_S2(:,bin)==fstr{bin}.pref_chain_S2(:,bin+6) & fstr{bin}.pref_chain_S2(:,bin)>0);
                end
            end
        end
%         save('0831_conn_mat_6s.mat','conn_mat_S1','conn_mat_S2');
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Generate data export for gephi network figure%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assume p_mat_sum(p_mat_sum<100)=0; c_mat_sum(p_mat_sum<100)=0; ckeep; in memory

% csvcell(1,:)={'Source','Target','Weight','Partition','Ratio'};
csvcell=cell(0,5);
for i=1:(length(pair_mat)-1)
    for j=(i+1):length(pair_mat)
        if sum(pair_mat(i,j,:))<50 %equivalent of 100 in two directions
%             csvcell(end+1,:)={reg_set{i},reg_set{j},1,1,0};
            continue
        else
            s1c=sum(conn_mat_S1(i,j,:))+sum(conn_mat_S1(j,i,:));
            s2c=sum(conn_mat_S2(i,j,:))+sum(conn_mat_S2(j,i,:));
            part=(s1c>0)*2+(s2c>0)*4;
            csvcell(end+1,:)={reg_set{i},reg_set{j},s1c+s2c,part,(s1c+s2c)/2/sum(pair_mat(i,j,:))};
        end
    end
end
writecell(csvcell,'selec_sum_ratio.csv');
    
end
