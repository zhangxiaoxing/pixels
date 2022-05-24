fstr=cell(1,6);
congru_all=[];

for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    t=fstr{bin}.conn_chain_S1(:,[1,2,2,2]);
    t(:,3)=1;
    t(:,4)=bin;
    t=[t,fstr{bin}.reg_chain_S1];
    sel=fstr{bin}.pref_chain_S1(:,bin)==fstr{bin}.pref_chain_S1(:,bin+6) & fstr{bin}.pref_chain_S1(:,bin)>0;
    congru_all=[congru_all;t(sel,:)];
    
    t=fstr{bin}.conn_chain_S2(:,[1,2,2,2]);
    t(:,3)=2;
    t(:,4)=bin;
    t=[t,fstr{bin}.reg_chain_S2];
    sel=fstr{bin}.pref_chain_S2(:,bin)==fstr{bin}.pref_chain_S2(:,bin+6) & fstr{bin}.pref_chain_S2(:,bin)>0;
    congru_all=[congru_all;t(sel,:)];

end

suids=unique([congru_all(:,1);congru_all(:,2)]);
corr_list=[];
for i=reshape(suids,1,[])
    su_reg=unique(congru_all(congru_all(:,1)==i,[1 2 6]),'rows');
    reg_su=unique(congru_all(congru_all(:,2)==i,[1 2 5]),'rows');
    
    regone=reshape(unique([su_reg(:,3);reg_su(:,3)]),1,[]);
    for i=regone
        corr_list=[corr_list;nnz(su_reg(:,3)==i),nnz(reg_su(:,3)==i)];
    end
    
end

fh=figure('Color','w','Position',[100,100,235,235]);
sh=scatter(corr_list(:,1),corr_list(:,2),14,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.15);
[r p]=corr(corr_list(:,1),corr_list(:,2));
legend(sh,sprintf('r = %.3f,p = %.3f',r,p));
xlim([0,30])
ylim([0,30])
set(gca,'XTick',0:10:30,'YTick',0:10:30)
xlabel('Memory neuron to region connection');
ylabel('Region to memory neuron connection');
exportgraphics(fh,'su_reg_recip_poo_plot.pdf')
