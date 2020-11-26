%% common dataset
fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
load('reg_keep.mat','reg_set')
sess=99;
lbound=sess*100000;
ubound=lbound+100000;

%% congruent
suid_all=[];
reg_all=[];
conn_all=[];
for bin=1:6
    sess_sel_s1=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound;
    congru_sel_s1=fstr{bin}.pref_chain_S1(:,bin)==1 & fstr{bin}.pref_chain_S1(:,bin+6)==1;
    sessconn=fstr{bin}.conn_chain_S1(sess_sel_s1 & congru_sel_s1,:);
%     [suids,BB,CC]=unique(sessconn);
    suids=sessconn(:);
    treg=fstr{bin}.reg_chain_S1(sess_sel_s1 & congru_sel_s1,:);
    reg=treg(:);
    suid_all=[suid_all;suids];
    reg_all=[reg_all;reg];
    conn_all=[conn_all;sessconn];
end
[suids,sidx,~]=unique(suid_all);
reg=reg_all(sidx);
[reg,I]=sort(reg);
suids=suids(I);
congru_conn_mat=zeros(length(suids),length(suids));
[GC,GR]=groupcounts(reg);
%     grpsel=reg<116 & ismember(reg,GR(GC>=15));
% [{reg_set{GR(GC>=15)}}',num2cell(GC(GC>=15))]
%     suids=suids(grpsel);
%     reg=reg(grpsel);

for i=1:length(suids)
    for j=1:length(suids)
        if i==j
            continue
        end
        if any(conn_all(:,1)==suids(i) & conn_all(:,2)==suids(j))
            if reg(i)==reg(j)
                congru_conn_mat(i,j)=1;
            else
                congru_conn_mat(i,j)=2;
            end
        end
    end
end


%%
cmap=[1,1,1;0.5,0.5,0.5;1,0,0];
fh=figure('Color','w','Position',[100,100,215,215]);
hold on
imagesc(congru_conn_mat)
colormap(cmap);

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
