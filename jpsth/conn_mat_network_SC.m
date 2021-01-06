%% common dataset

fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
load('reg_keep.mat','reg_set')

load reg_coord.mat
%% session loop
sess=99
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
%     congru_conn_mat_S2=zeros(length(suids_S2),length(suids_S2));
% 
%     for i=1:length(suids_S2)
%         for j=1:length(suids_S2)
%             if i==j
%                 continue
%             end
%             if any(conn_all_S2(:,1)==suids_S2(i) & conn_all_S2(:,2)==suids_S2(j))
%                 if reg_S2(i)~=reg_S2(j)
%                     congru_conn_mat_S2(i,j)=4;
%                 end
%             end
%         end
%     end
%     

    
    %% incongruent
    suid_all=[];
    reg_all=[];
    conn_all_incong=[];

    for bin=1:6
        sess_sel_s1=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound;
        incongru_sel_s1=diff(fstr{bin}.pref_chain_S1(:,[bin,bin+6]),1,2)~=0 & min(fstr{bin}.pref_chain_S1(:,[bin,bin+6]),[],2)>0;
        sessconn=fstr{bin}.conn_chain_S1(sess_sel_s1 & incongru_sel_s1,:);
        suids=sessconn(:);
        treg=fstr{bin}.reg_chain_S1(sess_sel_s1 & incongru_sel_s1,:);
        reg=treg(:);
        suid_all=[suid_all;suids];
        reg_all=[reg_all;reg];
        conn_all_incong=[conn_all_incong;sessconn];
        
        sess_sel_S2=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound;
        incongru_sel_S2=diff(fstr{bin}.pref_chain_S2(:,[bin,bin+6]),1,2)~=0 & min(fstr{bin}.pref_chain_S2(:,[bin,bin+6]),[],2)>0;
        sessconn=fstr{bin}.conn_chain_S2(sess_sel_S2 & incongru_sel_S2,:);
        suids=sessconn(:);
        treg=fstr{bin}.reg_chain_S2(sess_sel_S2 & incongru_sel_S2,:);
        reg=treg(:);
        suid_all=[suid_all;suids];
        reg_all=[reg_all;reg];
        conn_all_incong=[conn_all_incong;sessconn];
        
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
    
    suids_all=[suids_S1;suids_S2];
    sels_all=[zeros(size(suids_S1));ones(size(suids_S2))];
    reg_all=[reg_S1;reg_S2];
    [GC,GR]=groupcounts(reg_all);
    
    grp_sel=ismember(reg_all,GR(GC>15));
    
    reg_all=reg_all(grp_sel);
    sels_all=sels_all(grp_sel);
    suids_all=suids_all(grp_sel);
    
    [~,sidx]=sort(reg_all*1000+sels_all);
    suids_all=suids_all(sidx);
    reg_all=reg_all(sidx);
    
    conn_mat=zeros(length(suids_all),length(suids_all));
    
    for i=1:length(suids_all)
        for j=1:length(suids_all)
            if reg_all(i)==reg_all(j)
                continue
            end

            if any(conn_all_S1(:,1)==suids_all(i) & conn_all_S1(:,2)==suids_all(j))
                conn_mat(i,j)=1;
            end
            if any(conn_all_S2(:,1)==suids_all(i) & conn_all_S2(:,2)==suids_all(j))
                conn_mat(i,j)=2;
            end
            
            if any(conn_all_incong(:,1)==suids_all(i) & conn_all_incong(:,2)==suids_all(j))
                conn_mat(i,j)=3;
            end
            
        end
    end
    


    cmap=ones(4,3);
    cmap(2,:)=[1,0,0];
    cmap(3,:)=[1,0,0];
    cmap(4,:)=[0,0,1];

    

    fh=figure('Color','w','Position',[100,100,215,215]);
    hold on
    imagesc(conn_mat,[0,3]);
    colormap(cmap);






arrayfun(@(x) xline(x,':k'),[16,38,60,85,103]+0.5);
arrayfun(@(x) yline(x,':k'),[16,38,60,85,103]+0.5);
set(gca(),'XTick',[8.5,27,48,72.5,94,121],'XTickLabel',{'CA1','DP','ILA','MOs','PL','TTd'},...
    'YTick',[8.5,27,48,72.5,94,121],'YTickLabel',{'CA1','DP','ILA','MOs','PL','TTd'})
fh.Position(3:4)=[215,215];
xlim([0,143.5])
ylim([0,143.5])
exportgraphics(fh,'conn_mat_showcase_sess_99.pdf');

