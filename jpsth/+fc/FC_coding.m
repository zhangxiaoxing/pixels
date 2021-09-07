% assumes all stats file in memory for data generation
denovo=true;
memory=false;
keyboard();
%% Memory
if exist('denovo','var') && denovo
    su_var_t=[];
    keymap=containers.Map('KeyType','uint64','ValueType','uint64');
    for bin=1:6
        disp(bin);
        if memory
            load(sprintf('0116_memory_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        else
            load(sprintf('0115_nonsel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        end
        %% TODO: brain region
        for i=1:length(conn_chain_S1)
            if reg_chain_S1(i,1)==reg_chain_S1(i,2)
                continue
            end
            if rem(size(su_var_t,1),10000)==0
                disp(size(su_var_t,1))
            end
            key=conn_chain_S1(i,1)*100000+rem(conn_chain_S1(i,2),100000);
            if keymap.isKey(key)
                idx=keymap(key);
            else
                su_var_t(end+1,:)=zeros(1,13);
                su_var_t(end,1)=key;
                idx=size(su_var_t,1);
                keymap(key)=idx;
            end
            if memory
                if pref_chain_S1(i,bin) == pref_chain_S1(i,bin+6) && pref_chain_S1(i,bin)>0
                    selType=3;
                elseif max(pref_chain_S1(i,1:6)) == max(pref_chain_S1(i,7:12)) && max(pref_chain_S1(i,1:6))>0
                    selType=2;
                else
                    selType=1;
                end
            else
                selType=1;
            end
            su_var_t(idx,bin+1)=selType;
        end
        
        for i=1:length(conn_chain_S2)
            if reg_chain_S2(i,1)==reg_chain_S2(i,2)
                continue
            end
            if rem(size(su_var_t,1),10000)==0
                disp(size(su_var_t,1))
            end
            key=conn_chain_S2(i,1)*100000+rem(conn_chain_S2(i,2),100000);
            if keymap.isKey(key)
                idx=keymap(key);
            else
                su_var_t(end+1,:)=zeros(1,13);
                su_var_t(end,1)=key;
                idx=size(su_var_t,1);
                keymap(key)=idx;
            end
            if memory
                if pref_chain_S2(i,bin) == pref_chain_S2(i,bin+6) && pref_chain_S2(i,bin)>0
                    selType=3;
                elseif max(pref_chain_S2(i,1:6)) == max(pref_chain_S2(i,7:12)) && max(pref_chain_S2(i,1:6))>0
                    selType=2;
                else
                    selType=1;
                end
            else
                selType=1;
            end
            su_var_t(idx,bin+7)=selType;
        end
    end
    if memory
        save('su_var_t_S1_S2.mat','su_var_t');
    else
        save('su_var_t_S1_S2_nonsel.mat','su_var_t');
    end
else
    load('su_var_t_S1_S2.mat','su_var_t');
end


%% TODO: sort 
sort_mat=su_var_t(:,2:end)>0;
sel_mat=su_var_t(:,2:end)>1;
% sortsum=sort_mat*(2.^(11:-1:6)')+sel_mat*(2.^(5:-1:0)');
sortsum=sort_mat(:,1:6)*(2.^(0:5)')-sort_mat(:,7:12)*(2.^(0:5)');

[~,sumI]=sort(sortsum,'descend');
fh=figure('Color','w','Position',[32,32,500,1000]);
for subidx=1:2
subplot(1,2,subidx)
imagesc(su_var_t(sumI,(1:6)+1+(subidx-1)*6));
% cmap=[1 1 1;0 0 1;0 1 1;1 0 0];
cmap=[1 1 1;1 0 0;1 0 0;1 0 0];
colormap(cmap)
xlabel('Time (1-sec bin)')
ylabel('Coupled memory neuron pair #')
if memory
    set(gca(),'YTick',0:100000:300000,'YTickLabel',{'0','100K','200K','300K'},'XTick',1:6);
else
    set(gca(),'YTick',0:2e5:12e5,'YTickLabel',{'0','200K','400K','600K','800K','1M','1.2M'},'XTick',1:6);
end
ylim([-0.5,max(ylim())]);
end
if memory
    exportgraphics(fh,'FC_selective_S1_S2.pdf','ContentType','vector')
else
    exportgraphics(fh,'FC_selective_S1_S2_nonsel.pdf','ContentType','vector')
end




%% session showcase

%% common dataset

fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0116_memory_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
load('reg_keep.mat','reg_set')

load reg_coord.mat

sess=83;
lbound=sess*100000;
ubound=lbound+100000;
csvcell=cell(0,4);    
for bin=1:6
    sess_sel_S1=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound;
    reg_sel=fstr{bin}.reg_chain_S1(:,1)~=fstr{bin}.reg_chain_S1(:,2);
    sessconn=fstr{bin}.conn_chain_S1(sess_sel_S1 & reg_sel,:);
    for i=1:size(sessconn,1)
        csvcell(end+1,:)={sessconn(i,1),sessconn(i,2),bin,1};
    end
    
    sess_sel_S2=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound;
    reg_sel=fstr{bin}.reg_chain_S2(:,1)~=fstr{bin}.reg_chain_S2(:,2);
    sessconn=fstr{bin}.conn_chain_S2(sess_sel_S2 & reg_sel,:);
    for i=1:size(sessconn,1)
        csvcell(end+1,:)={sessconn(i,1),sessconn(i,2),bin,2};
    end
end

csvcell=[{'Source','Target','Bin','Partition'};csvcell];

writecell(csvcell,'FC_coding.csv');
