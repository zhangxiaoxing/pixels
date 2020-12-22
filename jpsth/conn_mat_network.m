keyboard();
if false
    congru_S1=readmatrix('congru_S1_conn_inter_sums.csv','Whitespace',' []');
    congru_S2=readmatrix('congru_S2_conn_inter_sums.csv','Whitespace',' []');
    incongru=readmatrix('incongru_conn_inter_sums.csv','Whitespace',' []');
    comb_S1=readmatrix('comb_S1_conn_conn_inter_sums.csv','Whitespace',' []');
    comb_S2=readmatrix('comb_S2_conn_conn_inter_sums.csv','Whitespace',' []');
    congru=[congru_S1;congru_S2];
    incongru=repmat(incongru,2,1);
    comb=[comb_S1;comb_S2];
    
    node_sel=congru(:,1)>20 & incongru(:,1)>20;
    one_comp=node_sel & congru(:,3)==1 & incongru(:,3)==1;
    
    plotOne(congru(node_sel,1),incongru(node_sel,1),'Node');
    plotOne(congru(node_sel,6),incongru(node_sel,6),'Avg. degree');
    plotOne(congru(node_sel,4),incongru(node_sel,4),'Avg. cluster coef.');
    
    plotOne(congru(node_sel,4),comb(node_sel,4),'Avg. cc. COMB');
    
    plotOne(congru(node_sel,5),incongru(node_sel,5),'Density');
    
    plotOne(congru(node_sel,7),incongru(node_sel,7),'Avg. path length');
    plotOne(congru(node_sel,3),incongru(node_sel,3),'Connected comp.');
    plotOne(congru(one_comp,7),incongru(one_comp,7),'1Comp path length');
end

%% common dataset

fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
load('reg_keep.mat','reg_set')

load reg_coord.mat
%% session loop
for sess=1:114
    lbound=sess*100000;
    ubound=lbound+100000;
    
    %% congruent S1
    suid_all_S1=[];
    reg_all_S1=[];
    conn_all_S1=[];
    for bin=1:6
        sess_sel_S1=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound;
        congru_sel_S1=fstr{bin}.pref_chain_S1(:,bin)==1 & fstr{bin}.pref_chain_S1(:,bin+6)==1;
        sessconn=fstr{bin}.conn_chain_S1(sess_sel_S1 & congru_sel_S1,:);
        suids=sessconn(:);
        treg=fstr{bin}.reg_chain_S1(sess_sel_S1 & congru_sel_S1,:);
        reg=treg(:);
        suid_all_S1=[suid_all_S1;suids];
        reg_all_S1=[reg_all_S1;reg];
        conn_all_S1=[conn_all_S1;sessconn];
    end
    
    [suids_S1,sidx_S1,~]=unique(suid_all_S1);
    reg_S1=reg_all_S1(sidx_S1);
    [reg_S1,I]=sort(reg_S1);
    suids_S1=suids_S1(I);
    grpsel_S1=reg_S1<116;% & ismember(reg,GR(GC>=15));

    suids_S1=suids_S1(grpsel_S1);
    reg_S1=reg_S1(grpsel_S1);
    [GC_S1,GR_S1]=groupcounts(reg_S1);
    congru_conn_mat_S1=zeros(length(suids_S1),length(suids_S1));
%     congru_conn_mat_S2=zeros(length(suids_S2),length(suids_S2));
    for i=1:length(suids_S1)
        for j=1:length(suids_S1)
            if i==j
                continue
            end
            if any(conn_all_S1(:,1)==suids_S1(i) & conn_all_S1(:,2)==suids_S1(j))
                if reg_S1(i)~=reg_S1(j)
                    congru_conn_mat_S1(i,j)=4;
                end
            end
        end
    end
    
    congru_mat_S1=congru_conn_mat_S1;
    csvcell_S1={'Source','Target','Pref'};
    for i=1:length(congru_mat_S1)
        for j=1:length(congru_mat_S1)
            if i==j || congru_mat_S1(i,j)==0
                continue
            end
            if congru_mat_S1(i,j)>=4
                csvcell_S1(end+1,:)={suids_S1(i),suids_S1(j),congru_mat_S1(i,j)};
            end
        end
    end
    writecell(csvcell_S1,sprintf('congru_conn_sc_inter_%03d_S1.csv',sess))
    
    csvcell_S1={'Id','Label','Reg','AP','DV'};
    for i=1:length(congru_mat_S1)
        regidx_S1=cellfun(@(x) strcmp(x,reg_set{reg_S1(i)}),regcoord(:,1));
        csvcell_S1(end+1,:)={suids_S1(i),reg_set{reg_S1(i)},reg_set{reg_S1(i)},regcoord{regidx_S1,2}/15+rand(1)*2-1,1000/15-regcoord{regidx_S1,3}/15+rand(1)*2-1};
    end
    writecell(csvcell_S1,sprintf('congru_conn_sc_node_coord_%03d_S1.csv',sess))
    
    %% congruent S2
    suid_all_S2=[];
    reg_all_S2=[];
    conn_all_S2=[];
    for bin=1:6
        sess_sel_S2=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound;
        congru_sel_S2=fstr{bin}.pref_chain_S2(:,bin)==2 & fstr{bin}.pref_chain_S2(:,bin+6)==2;
        sessconn=fstr{bin}.conn_chain_S2(sess_sel_S2 & congru_sel_S2,:);
        suids=sessconn(:);
        treg=fstr{bin}.reg_chain_S2(sess_sel_S2 & congru_sel_S2,:);
        reg=treg(:);
        suid_all_S2=[suid_all_S2;suids];
        reg_all_S2=[reg_all_S2;reg];
        conn_all_S2=[conn_all_S2;sessconn];
    end
    
    [suids_S2,sidx_S2,~]=unique(suid_all_S2);
    reg_S2=reg_all_S2(sidx_S2);
    [reg_S2,I]=sort(reg_S2);
    suids_S2=suids_S2(I);
    grpsel_S2=reg_S2<116;% & ismember(reg,GR(GC>=15));

    suids_S2=suids_S2(grpsel_S2);
    reg_S2=reg_S2(grpsel_S2);
    [GC_S2,GR_S2]=groupcounts(reg_S2);
    congru_conn_mat_S2=zeros(length(suids_S2),length(suids_S2));
%     congru_conn_mat_S2=zeros(length(suids_S2),length(suids_S2));
    for i=1:length(suids_S2)
        for j=1:length(suids_S2)
            if i==j
                continue
            end
            if any(conn_all_S2(:,1)==suids_S2(i) & conn_all_S2(:,2)==suids_S2(j))
                if reg_S2(i)~=reg_S2(j)
                    congru_conn_mat_S2(i,j)=4;
                end
            end
        end
    end
    
    congru_mat_S2=congru_conn_mat_S2;
    csvcell_S2={'Source','Target','Pref'};
    for i=1:length(congru_mat_S2)
        for j=1:length(congru_mat_S2)
            if i==j || congru_mat_S2(i,j)==0
                continue
            end
            if congru_mat_S2(i,j)>=4
                csvcell_S2(end+1,:)={suids_S2(i),suids_S2(j),congru_mat_S2(i,j)};
            end
        end
    end
    writecell(csvcell_S2,sprintf('congru_conn_sc_inter_%03d_S2.csv',sess))
    
    csvcell_S2={'Id','Label','Reg','AP','DV'};
    for i=1:length(congru_mat_S2)
        regidx_S2=cellfun(@(x) strcmp(x,reg_set{reg_S2(i)}),regcoord(:,1));
        csvcell_S2(end+1,:)={suids_S2(i),reg_set{reg_S2(i)},reg_set{reg_S2(i)},regcoord{regidx_S2,2}/15+rand(1)*2-1,1000/15-regcoord{regidx_S2,3}/15+rand(1)*2-1};
    end
    writecell(csvcell_S2,sprintf('congru_conn_sc_node_coord_%03d_S2.csv',sess))
    
    %% incongruent
    suid_all=[];
    reg_all=[];
    conn_all_S1=[];
    conn_all_S2=[];
    for bin=1:6
        sess_sel_s1=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound;
        incongru_sel_s1=diff(fstr{bin}.pref_chain_S1(:,[bin,bin+6]),1,2)~=0 & min(fstr{bin}.pref_chain_S1(:,[bin,bin+6]),[],2)>0;
        sessconn=fstr{bin}.conn_chain_S1(sess_sel_s1 & incongru_sel_s1,:);
        suids=sessconn(:);
        treg=fstr{bin}.reg_chain_S1(sess_sel_s1 & incongru_sel_s1,:);
        reg=treg(:);
        suid_all=[suid_all;suids];
        reg_all=[reg_all;reg];
        conn_all_S1=[conn_all_S1;sessconn];
        
        sess_sel_S2=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound;
        incongru_sel_S2=diff(fstr{bin}.pref_chain_S2(:,[bin,bin+6]),1,2)~=0 & min(fstr{bin}.pref_chain_S2(:,[bin,bin+6]),[],2)>0;
        sessconn=fstr{bin}.conn_chain_S2(sess_sel_S2 & incongru_sel_S2,:);
        suids=sessconn(:);
        treg=fstr{bin}.reg_chain_S2(sess_sel_S2 & incongru_sel_S2,:);
        reg=treg(:);
        suid_all=[suid_all;suids];
        reg_all=[reg_all;reg];
        conn_all_S2=[conn_all_S2;sessconn];
        
    end
    [suids,sidx,~]=unique(suid_all);
    reg=reg_all(sidx);
    [reg,I]=sort(reg);
    suids=suids(I);
    [GC,GR]=groupcounts(reg);
    grpsel=reg<116;% & ismember(reg,GR(GC>=15));
    % [{reg_set{GR(GC>=15)}}',num2cell(GC(GC>=15))]
    suids=suids(grpsel);
    reg=reg(grpsel);
    
    incongru_conn_mat_S1=zeros(length(suids),length(suids));
    incongru_conn_mat_S2=zeros(length(suids),length(suids));
    for i=1:length(suids)
        for j=1:length(suids)
            if i==j
                continue
            end
            if any(conn_all_S1(:,1)==suids(i) & conn_all_S1(:,2)==suids(j))
                if reg(i)==reg(j)
                    incongru_conn_mat_S1(i,j)=1;
                else
                    incongru_conn_mat_S1(i,j)=4;
                end
            end
            
            if any(conn_all_S2(:,1)==suids(i) & conn_all_S2(:,2)==suids(j))
                if reg(i)==reg(j)
                    incongru_conn_mat_S2(i,j)=2;
                else
                    incongru_conn_mat_S2(i,j)=8;
                end
            end
            
        end
    end
    
    
    %%
%     cmap=ones(15,3);
%     cmap([1,4,5],:)=repmat([0.5],3,3);
%     cmap(2,:)=[0.8,0,0];
%     cmap(8,:)=[0,0,0.8];
%     cmap(10,:)=[1,0,1];
%     cmap=[1,1,1;cmap];
    incongru_mat=incongru_conn_mat_S1+incongru_conn_mat_S2;
%     fh=figure('Color','w','Position',[100,100,215,215]);
%     hold on
%     imagesc(incongru_mat,[0,15]);
%     colormap(cmap);
%     xlim([0,length(congru_conn_mat_S1)])
%     ylim([0,length(congru_conn_mat_S1)])
    csvcell={'Source','Target','Pref'};
    for i=1:length(incongru_mat)
        for j=1:length(incongru_mat)
            if i==j || incongru_mat(i,j)==0
                continue
            end
            if incongru_mat(i,j)<4
%                 csvcell(end+1,:)={i,j,1};
            else
                csvcell(end+1,:)={i,j,12};
            end
        end
    end
    writecell(csvcell,sprintf('incongru_conn_sc_inter_%03d.csv',sess))
    
    csvcell={'Id','Label','Reg','AP','DV'};
    for i=1:length(incongru_mat)
        regidx=cellfun(@(x) strcmp(x,reg_set{reg(i)}),regcoord(:,1));
        csvcell(end+1,:)={i,reg_set{reg(i)},reg_set{reg(i)},regcoord{regidx,2}/15+rand(1)*2-1,1000/15-regcoord{regidx,3}/15+rand(1)*2-1};
    end
    writecell(csvcell,sprintf('incongru_conn_sc_node_coord_%03d.csv',sess))
    
end
return


















arrayfun(@(x) xline(x,':k'),[16,38,59,85,102,139]+0.5);
arrayfun(@(x) yline(x,':k'),[16,38,59,85,102,139]+0.5);
set(gca(),'XTick',[8.5,27,48,72.5,94,121],'XTickLabel',{'CA1','DP','ILA','MOs','PL','TTd'},...
    'YTick',[8.5,27,48,72.5,94,121],'YTickLabel',{'CA1','DP','ILA','MOs','PL','TTd'})
fh.Position(3:4)=[215,215];
xlim([0,139.5])
ylim([0,139.5])
exportgraphics(fh,'conn_mat_showcase_sess_99.pdf');

pair_sess_sel=fstr.pair_chain(:,1)>=lbound & fstr.pair_chain(:,1)<ubound;
sess_pair_reg=fstr.pair_reg(pair_sess_sel,:);

tsu=fstr.pair_chain(pair_sess_sel,:);
reglist=unique(reg);
% member_pair_reg=sess_pair_reg(all(ismember(sess_pair_reg,reglist),2),:);
% numel(unique(tsu(ismember(sess_pair_reg,reglist))));
for i=1:numel(reglist)
    for j=1:numel(reglist)
        if i==j
            continue
        else
            conn_count=nnz(treg(:,1)==reglist(i) & treg(:,2)==reglist(j));
            pair_count=nnz((sess_pair_reg(:,1)==reglist(i) & sess_pair_reg(:,2)==reglist(j)) ...
                | (sess_pair_reg(:,1)==reglist(j) & sess_pair_reg(:,2)==reglist(i)));
            
            disp({reg_set{reglist(i)},reg_set{reglist(j)},conn_count,pair_count})
        end
    end
end



return


inter_reg_cong_count=zeros(114,1);
sess_reg=cell(114,1);
for i=1:114
    lbound=100000*i;
    ubound=100000*(i+1);
    sel=fstr{6}.conn_chain_S1(:,1)>lbound & fstr{6}.conn_chain_S1(:,1)<ubound;
    inter_reg_cong_count(i)=nnz(diff(fstr{6}.reg_chain_S1(sel,:),1,2) & fstr{6}.pref_chain_S1(sel,6)==fstr{6}.pref_chain_S1(sel,12));
    sess_reg{i}=arrayfun(@(x) reg_set{x}, unique(fstr{6}.reg_chain_S1(sel,:)),'UniformOutput',false);
end
[C,I]=sort(inter_reg_cong_count,'descend');
sess_reg=[C,sess_reg(I)];


function plotOne(congru,incongru,lbl)
fh=figure('Color','w','Position',[100,100,125,200]);
hold on;
ydata=[incongru,congru];
xdata=ones(size(ydata));
if contains(lbl,'comp.')
    xdata(:,1)=xdata(:,1)+rand(size(ydata,1),1)*0.3;
    xdata(:,2)=xdata(:,2)+1-rand(size(ydata,1),1)*0.3;
else
    xdata(:,1)=xdata(:,1)+rand(size(ydata,1),1)*0.2;
    xdata(:,2)=xdata(:,2)+1-rand(size(ydata,1),1)*0.2;
end
plot(xdata',ydata','-k')
plot(xdata(:,1),ydata(:,1),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','none')
plot(xdata(:,2),ydata(:,2),'ko','MarkerSize',4,'MarkerFaceColor','none','MarkerEdgeColor','k')
errorbar([0.75,2.25],mean(ydata),std(ydata)/sqrt(size(ydata,1)),'ko','MarkerSize',6);
xlim([0.5,2.5])
ylim([0,max(ylim())])
ylabel(lbl)
xlabel(sprintf('%0.4f',signrank(congru,incongru)));
set(gca,'XTick',[])
if contains(lbl,'comp.')
    set(gca,'YScale','log')
    ylim([0.9,11])
    ax=gca;
    ax.YAxis.MinorTickValues=1:9;
    ax.YAxis.MinorTick='on';
    ax.TickLength=[0.03,0.06]
end

exportgraphics(fh,sprintf('func_conn_stats_%s.pdf',lbl),'ContentType','vector');
end
