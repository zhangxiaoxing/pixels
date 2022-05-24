clear

%% Prepartion
if false
    load('conn_nonsel')
    load('conn.mat')
    load('D:\code\FR_modulated_20210104.mat','FR_modulated')
    load('D:\code\114_sorted_file_path.mat')
    
    [conn_congru_active,list_congru_active]=PlotList(congru_active_all,FR_modulated,sorted_fpath);
    [congru_active_pair,list_congru_active_pair]=PlotList(congru_active_pair_all,FR_modulated,sorted_fpath);
    [conn_congru_inactive,list_congru_inactive]=PlotList(congru_inactive_all,FR_modulated,sorted_fpath);
    [congru_inactive_pair,list_congru_inactive_pair]=PlotList(congru_inactive_pair_all,FR_modulated,sorted_fpath);
    [conn_incong_active,list_incong_active]=PlotList(incong_active_all,FR_modulated,sorted_fpath);
    [incong_active_pair,list_incong_active_pair]=PlotList(incong_active_pair_all,FR_modulated,sorted_fpath);
    [conn_incong_inactive,list_incong_inactive]=PlotList(incong_inactive_all,FR_modulated,sorted_fpath);
    [incong_inactive_pair,list_incong_inactive_pair]=PlotList(incong_inactive_pair_all,FR_modulated,sorted_fpath);
    [conn_nonsel,list_nonsel]=PlotList(conn_nonse,FR_modulated,sorted_fpath);
    [nonsel_pair,list_nonsel_pair]=PlotList(non_pair,FR_modulated,sorted_fpath);
    
    save('D:/code/conn_FR_list_20200104.mat','list_congru_active','list_congru_inactive','list_incong_active','list_incong_inactive','list_nonsel','conn_congru_active','conn_congru_inactive','conn_incong_active','conn_incong_inactive','conn_nonsel','-v7.3')
    save('D:/code/pair_FR_list_20200104.mat','list_congru_active_pair','list_congru_inactive_pair','list_incong_active_pair','list_incong_inactive_pair','list_nonsel_pair','congru_active_pair','congru_inactive_pair','incong_active_pair','incong_inactive_pair','nonsel_pair','-v7.3')
end
%% Coupling ratio in different FR
if false
load('D:/code/pair_FR_list_20200104.mat')
load('D:/code/conn_FR_list_20200104.mat')
FR_all=[2,5,10];
e=0.5;
for f=1:length(FR_all)
    FR=FR_all(f);  
    
    r{f,5}=CalculatedRatio(list_congru_active,conn_congru_active,list_congru_active_pair,congru_active_pair,FR,e);
    r{f,3}=CalculatedRatio(list_incong_active,conn_incong_active,list_incong_active_pair,incong_active_pair,FR,e);
    r{f,4}=CalculatedRatio(list_congru_inactive,conn_congru_inactive,list_congru_inactive_pair,congru_inactive_pair,FR,e);
    r{f,2}=CalculatedRatio(list_incong_inactive,conn_incong_inactive,list_incong_inactive_pair,incong_inactive_pair,FR,e);
    r{f,1}=CalculatedRatio(list_nonsel,conn_nonsel,list_nonsel_pair,nonsel_pair,FR,e);

%     x=[connidx_congru_active(:,2);connidx_incong_active(:,2);connidx_congru_inactive(:,2);connidx_incong_inactive(:,2);connidx_nonsel(:,2)];
%     y=[ones(size(connidx_congru_active,1),1);repmat(2,size(connidx_incong_active,1),1);repmat(3,size(connidx_congru_inactive,1),1);repmat(4,size(connidx_incong_inactive,1),1);repmat(5,size(connidx_nonsel,1),1)];
%     p(f,1)=anova1(x,y);
    
    
end
close all
save('D:\code\conn_diffFR_0104.mat','r')
c={'k','r','b'};
for i=1:size(m,1)
    for j=1:5
    m(i,j)=mean(r{i,j});
    s(i,j)=std(r{i,j})/sqrt(size(r{i,j},1));
    end
   
end
 hold on
 for i=1:size(m,1)
     errorbar(1:5,m(i,:),s(i,:),'.-','CapSize',12,'Color',c{i})
 end
set(gca,'XTick',1:5,'XTickLabel',{'non-memory','incongruent inactive','incongruent activate','congruent inactive','congruent active'},'XTickLabelRotation',30,'FontSize',12,'FontName','Arial')
box off
xlim([0 6])
exportgraphics(gcf,sprintf('Conn_ratio_FR%d.pdf',0104),'ContentType','vector');

x=[];y1=[];y2=[];
for i=1:5
    for j=1:3
       x=[x;r{j,i}];
       y1=[y1;ones(size(r{j,i},1),1)*i];
       y2=[y2;ones(size(r{j,i},1),1)*j];
    end
end
anovan(x,{y1,y2})
end
%% Coupling ratio in FR-increased/decreased neurons
clear
load('D:/code/pair_FR_list_20200104.mat')
load('D:/code/conn_FR_list_20200104.mat')
%
conn_congru_active=FRchange(list_congru_active,conn_congru_active);
congru_active_pair=FRchange(list_congru_active_pair,congru_active_pair);
ratio=[];
for j=1:114
    ratio=[ratio;[nnz(conn_congru_active(:,1)>=j*100000&conn_congru_active(:,1)<(j+1)*100000 &conn_congru_active(:,4)==-1 &conn_congru_active(:,5)==-1)/2/nnz(congru_active_pair(:,1)>=j*100000&congru_active_pair(:,1)<(j+1)*100000 &congru_active_pair(:,4)==-1&congru_active_pair(:,5)==-1),...
        nnz(conn_congru_active(:,1)>=j*100000&conn_congru_active(:,1)<(j+1)*100000 &conn_congru_active(:,4)==1&conn_congru_active(:,5)==1)/2/nnz(congru_active_pair(:,1)>=j*100000&congru_active_pair(:,1)<(j+1)*100000 &congru_active_pair(:,4)==1&congru_active_pair(:,5)==1),...
        nnz(conn_congru_active(:,1)>=j*100000&conn_congru_active(:,1)<(j+1)*100000 &(conn_congru_active(:,4).*conn_congru_active(:,5)==-1))/2/nnz(congru_active_pair(:,1)>=j*100000&congru_active_pair(:,1)<(j+1)*100000 & (congru_active_pair(:,4).*congru_active_pair(:,5)==-1))]];
end
ratio(ratio(:,1)==Inf|ratio(:,2)==Inf|ratio(:,3)==Inf|isnan(ratio(:,1))|isnan(ratio(:,2))|isnan(ratio(:,3)),:)=[];
r{1,5}=ratio(ratio(:,1)~=0,1);
r{2,5}=ratio(ratio(:,2)~=0,2);
r{3,5}=ratio(ratio(:,3)~=0,3);

conn_incong_active=FRchange(list_incong_active,conn_incong_active);
incong_active_pair=FRchange(list_incong_active_pair,incong_active_pair);
ratio=[];
for j=1:114
    ratio=[ratio;[nnz(conn_incong_active(:,1)>=j*100000&conn_incong_active(:,1)<(j+1)*100000 &conn_incong_active(:,4)==-1&conn_incong_active(:,5)==-1)/2/nnz(incong_active_pair(:,1)>=j*100000&incong_active_pair(:,1)<(j+1)*100000 &incong_active_pair(:,4)==-1&incong_active_pair(:,5)==-1),...
        nnz(conn_incong_active(:,1)>=j*100000&conn_incong_active(:,1)<(j+1)*100000 &conn_incong_active(:,4)==1 &conn_incong_active(:,5)==1)/2/nnz(incong_active_pair(:,1)>=j*100000&incong_active_pair(:,1)<(j+1)*100000 &incong_active_pair(:,4)==1&incong_active_pair(:,5)==1),...
        nnz(conn_incong_active(:,1)>=j*100000&conn_incong_active(:,1)<(j+1)*100000 &conn_incong_active(:,4).*conn_incong_active(:,5)==-1)/2/nnz(incong_active_pair(:,1)>=j*100000&incong_active_pair(:,1)<(j+1)*100000 &incong_active_pair(:,4).*incong_active_pair(:,5)==-1)]];
end
ratio(ratio(:,1)==Inf|ratio(:,2)==Inf|ratio(:,3)==Inf|isnan(ratio(:,1))|isnan(ratio(:,2))|isnan(ratio(:,3)),:)=[];
r{1,3}=ratio(ratio(:,1)~=0,1);
r{2,3}=ratio(ratio(:,2)~=0,2);
r{3,3}=ratio(ratio(:,3)~=0,3);

conn_congru_inactive=FRchange(list_congru_inactive,conn_congru_inactive);
congru_inactive_pair=FRchange(list_congru_inactive_pair,congru_inactive_pair);
ratio=[];
for j=1:114
    ratio=[ratio;[nnz(conn_congru_inactive(:,1)>=j*100000&conn_congru_inactive(:,1)<(j+1)*100000 &conn_congru_inactive(:,4)==-1 &conn_congru_inactive(:,5)==-1)/2/nnz(congru_inactive_pair(:,1)>=j*100000&congru_inactive_pair(:,1)<(j+1)*100000 &congru_inactive_pair(:,4)==-1 &congru_inactive_pair(:,5)==-1),...
        nnz(conn_congru_inactive(:,1)>=j*100000&conn_congru_inactive(:,1)<(j+1)*100000 &conn_congru_inactive(:,4)==1 &conn_congru_inactive(:,5)==1)/2/nnz(congru_inactive_pair(:,1)>=j*100000&congru_inactive_pair(:,1)<(j+1)*100000 &congru_inactive_pair(:,4)==1 &congru_inactive_pair(:,5)==1),...
        nnz(conn_congru_inactive(:,1)>=j*100000&conn_congru_inactive(:,1)<(j+1)*100000 &conn_congru_inactive(:,4).*conn_congru_inactive(:,5)==-1)/2/nnz(congru_inactive_pair(:,1)>=j*100000&congru_inactive_pair(:,1)<(j+1)*100000 &congru_inactive_pair(:,4).*congru_inactive_pair(:,5)==-1)]];
end
ratio(ratio(:,1)==Inf|ratio(:,2)==Inf|ratio(:,3)==Inf|isnan(ratio(:,1))|isnan(ratio(:,2))|isnan(ratio(:,3)),:)=[];
r{1,4}=ratio(ratio(:,1)~=0,1);
r{2,4}=ratio(ratio(:,2)~=0,2);
r{3,4}=ratio(ratio(:,3)~=0,3);

conn_incong_inactive=FRchange(list_incong_inactive,conn_incong_inactive);
incong_inactive_pair=FRchange(list_incong_inactive_pair,incong_inactive_pair);ratio=[];
for j=1:114
    ratio=[ratio;[nnz(conn_incong_inactive(:,1)>=j*100000&conn_incong_inactive(:,1)<(j+1)*100000 &(conn_incong_inactive(:,4)==-1 & conn_incong_inactive(:,5)==-1))/2/nnz(incong_inactive_pair(:,1)>=j*100000&incong_inactive_pair(:,1)<(j+1)*100000 &incong_inactive_pair(:,4)==-1&incong_inactive_pair(:,5)==-1),...
        nnz(conn_incong_inactive(:,1)>=j*100000&conn_incong_inactive(:,1)<(j+1)*100000 &conn_incong_inactive(:,4)==1& conn_incong_inactive(:,5)== 1)/2/nnz(incong_inactive_pair(:,1)>=j*100000&incong_inactive_pair(:,1)<(j+1)*100000 &incong_inactive_pair(:,4)==1&incong_inactive_pair(:,5)==1),...
        nnz(conn_incong_inactive(:,1)>=j*100000&conn_incong_inactive(:,1)<(j+1)*100000 &conn_incong_inactive(:,4).* conn_incong_inactive(:,5)== -1)/2/nnz(incong_inactive_pair(:,1)>=j*100000&incong_inactive_pair(:,1)<(j+1)*100000 &incong_inactive_pair(:,4).*incong_inactive_pair(:,5)==-1)]];
end
ratio(ratio(:,1)==Inf|ratio(:,2)==Inf|ratio(:,3)==Inf|isnan(ratio(:,1))|isnan(ratio(:,2))|isnan(ratio(:,3)),:)=[];
r{1,2}=ratio(ratio(:,1)~=0,1);
r{2,2}=ratio(ratio(:,2)~=0,2);
r{3,2}=ratio(ratio(:,3)~=0,3);
%
conn_nonsel=FRchange(list_nonsel,conn_nonsel);
nonsel_pair=FRchange(list_nonsel_pair,nonsel_pair);
ratio=[];
for j=1:114
    ratio=[ratio;[nnz(conn_nonsel(:,1)>=j*100000&conn_nonsel(:,1)<(j+1)*100000 &conn_nonsel(:,4)==-1&conn_nonsel(:,5)==-1)/2/nnz(nonsel_pair(:,1)>=j*100000&nonsel_pair(:,1)<(j+1)*100000 &nonsel_pair(:,4)==-1&nonsel_pair(:,5)==-1),...
        nnz(conn_nonsel(:,1)>=j*100000&conn_nonsel(:,1)<(j+1)*100000 &conn_nonsel(:,4)==1&conn_nonsel(:,5)==1)/2/nnz(nonsel_pair(:,1)>=j*100000&nonsel_pair(:,1)<(j+1)*100000 &nonsel_pair(:,4)==1 &nonsel_pair(:,5)==1),...
        nnz(conn_nonsel(:,1)>=j*100000&conn_nonsel(:,1)<(j+1)*100000 &conn_nonsel(:,4).*conn_nonsel(:,5)==-1)/2/nnz(nonsel_pair(:,1)>=j*100000&nonsel_pair(:,1)<(j+1)*100000 &nonsel_pair(:,4).*nonsel_pair(:,5)==-1)]];
end
ratio(ratio(:,1)==Inf|ratio(:,2)==Inf|ratio(:,3)==Inf|isnan(ratio(:,1))|isnan(ratio(:,2))|isnan(ratio(:,3)),:)=[];
r{1,1}=ratio(ratio(:,1)~=0,1);
r{2,1}=ratio(ratio(:,2)~=0,2);
r{3,1}=ratio(ratio(:,3)~=0,3);

save('Conn_FRchange.mat','r')

for i=1:5
    m(1,i)=mean(r{1,i},1);
    m(2,i)=mean(r{2,i},1);
    m(3,i)=mean(r{3,i},1);
    s(1,i)=std(r{1,i})/sqrt(114);
    s(2,i)=std(r{2,i})/sqrt(114);
    s(3,i)=std(r{3,i})/sqrt(114);
end
bar(0.6:2:10,m(1,:),0.2,'w','Edgecolor','b')
hold on
bar(1.0:2:10,m(3,:),0.2,'w','Edgecolor','k')
hold on
bar(1.4:2:10,m(2,:),0.2,'w','Edgecolor','r')
hold on
%     legend('Decreased-Decreased','Increased-Increased','Increased-Decreased/Decreased-Increased''Location','northwest')
errorbar(0.6:2:10,m(1,:),s(1,:),'b.','CapSize',12)
hold on
errorbar(1.0:2:10,m(3,:),s(3,:),'k.','CapSize',12)
hold on
errorbar(1.4:2:10,m(2,:),s(2,:),'r.','CapSize',12)
set(gca,'XTick',1:2:10,'XTickLabel',{'non-memory','incongruent inactive','incongruent activate','congruent inactive','congruent active'},'XTickLabelRotation',30,'FontSize',12,'FontName','Arial')
ylabel('Connection density','FontSize',18,'FontName','Arial');
box off
exportgraphics(gcf,'Conn_ratio_FR_0104.pdf','ContentType','vector');

x=[];y1=[];y2=[];
for i=1:5
    for j=1:3
       x=[x;r{j,i}];
       y1=[y1;ones(size(r{j,i},1),1)*i];
       y2=[y2;ones(size(r{j,i},1),1)*j];
    end
end
anovan(x,{y1,y2})

%% Histgram of FR in memory neuons
clear
load('D:/code/conn_FR_list_20200104.mat')
fh=figure('Color','w','Position',[100,100,400,400]);
subplot(2,2,1)
[counts,centers]=hist(list_congru_active(:,2)-list_congru_active(:,3),20);
bar(centers, counts / sum(counts))
title('congru active')
box off
subplot(2,2,2)
[counts,centers]=hist(list_incong_active(:,2)-list_incong_active(:,3),20);
bar(centers, counts / sum(counts))
title('incongru active')
box off
subplot(2,2,3)
[counts,centers]=hist(list_congru_inactive(:,2)-list_congru_inactive(:,3),20);
bar(centers, counts / sum(counts))
title('congru inactive')
box off
subplot(2,2,4)
[counts,centers]=hist(list_incong_inactive(:,2)-list_incong_inactive(:,3),20);
bar(centers, counts / sum(counts))
title('incongru inactive')
box off
exportgraphics(fh,'hist_FR_0104.pdf','ContentType','vector');   

%% FR in active/inactive/no-memory neurons
clear
load('D:\code\114_sorted_file_path.mat')
load('D:\code\FR_modulated_20210104.mat')
sus_trans=h5read('D:\code\transient_6_0104.hdf5','/sus_trans');
cid_list=h5read('D:\code\transient_6_0104.hdf5','/cluster_id');
path_list=h5read('D:\code\transient_6_0104.hdf5','/path');
reg_list=h5read('D:\code\transient_6_0104.hdf5','/reg');
load('reg_keep.mat')
reg_list=regexp(reg_list,'(\w|\\|-)*','match','once');
reg_logic=ismember(deblank(reg_list),reg_set(1:115));

sel=sus_trans((sus_trans(:,1)==1|sus_trans(:,2)==1)&reg_logic,:);
FR_sel=FR_modulated((sus_trans(:,1)==1|sus_trans(:,2)==1)&reg_logic,1);
nonsel=sus_trans((sus_trans(:,1)~=1&sus_trans(:,2)~=1)&reg_logic,:);
FR_nonsel=FR_modulated((sus_trans(:,1)~=1&sus_trans(:,2)~=1)&reg_logic,1);
for i=1:size(sel,1)
    for bin=1:6
        FR(bin)=mean(mean(FR_sel{i}(:,4*bin+13:4*bin+16),2),1);
    end
  FR_active(i,1)=mean(FR(sel(i,8:13)>0));
  FR_inactive(i,1)=mean(FR(sel(i,8:13)==0));  
end
FR_inactive(isnan(FR_inactive))=[];
for i=1:size(nonsel,1)
    for bin=1:6
        FR(bin)=mean(mean(FR_nonsel{i}(:,4*bin+13:4*bin+16),2),1);
    end
 FR_nonsel_all(i,1)=mean(FR);
end

fh=figure('Color','w','Position',[100,100,150,300]);
boxplot([FR_nonsel_all;FR_inactive;FR_active],[ones(size(FR_nonsel_all,1),1)*1;ones(size(FR_inactive,1),1)*2;ones(size(FR_active,1),1)*3])   
set(gca,'XTick',1:3,'XTickLabel',{'non-memory','inactive','active'},'XTickLabelRotation',30,'FontSize',12,'FontName','Arial')
box off
xlim([0 4])
ylabel('Firing rate (Hz)')
exportgraphics(fh,'FR_activebin_0105.pdf','ContentType','vector');

anova1([FR_nonsel_all;FR_inactive;FR_active],[ones(size(FR_nonsel_all,1),1)*1;ones(size(FR_inactive,1),1)*2;ones(size(FR_active,1),1)*3])   


%% function
function [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list,homedir)
metaFolder=replace(folder,'\','/');
metaFolder=fullfile(fullfile(homedir,'DataSum'),metaFolder);
if isfolder(metaFolder)
    spkFolder=replace(metaFolder,'imec1','imec0');
    file=dir(fullfile(spkFolder,'spike_info.mat'));
    if isempty(file)
        folderType=-1;
        file=[];
        spkFolder=[];
        disp('Error processing file 2-tracks');
        disp(metaFolder);
        error_list(end+1,:)={folderType,metaFolder};
        %             pause;
        return
    end
    folderType=2;
else
    metaFolder=replace(metaFolder,'DataSum','DataSum/singleProbe');
    spkFolder=metaFolder;
    file=dir(fullfile(spkFolder,'spike_times.npy'));
    if isempty(file)
        folderType=-1;
        file=[];
        spkFolder=[];
        disp('Error processing file 1-track');
        disp(metaFolder);
        error_list(end+1,:)={folderType,metaFolder};
        return
    end
    folderType=1;
end
end

function [conn,out]=PlotList(conn0,FR_modulated,sorted_fpath)
cid_list=double(h5read('D:\code\transient_6_0104.hdf5','/cluster_id'));
path_list=h5read('D:\code\transient_6_0104.hdf5','/path');
path_list=regexp(path_list,'(\w|\\|-)*','match','once');
reg_list=h5read('D:\code\transient_6_0104.hdf5','/reg');
load('D:\code\reg_keep.mat')
reg_list=regexp(reg_list,'(\w|\\|-)*','match','once');
reg_logic=ismember(deblank(reg_list),reg_set(1:115));

conn=[];

for bin=1:6
    for j=1:size(conn0,1)
        for i=1:size(conn0{j,bin},1)
            conn=[conn;[conn0{j,bin}{i,1},repmat(bin,size(conn0{j,bin}{i,1},1),1)]];
        end
    end
end
out(:,1)=unique(conn(:,1:2));

for i=1:size(out(:,1),1)
    sessIdx=round(out(i,1)/100000);
    folder=sorted_fpath{sessIdx};    
    [folderType,~,~,~,~]=jointFolder(folder,cell(0),'I:\WT');
    suidx=out(i,1)-sessIdx*100000;
    if folderType>1 && suidx>=10000
        folder=replace(folder,'imec0','imec1');
        suidx=suidx-10000;
    elseif folderType>1 && suidx<=10000
        folder=replace(folder,'imec1','imec0');
    end
    try 
    FR=FR_modulated{contains(path_list,folder)&cid_list==suidx&reg_logic,1};
    catch 
        continue
    end   
    out(i,2)=mean(mean(FR(17:40),2),1);
    out(i,3)=mean(mean(FR(1:10),2),1); 
end

out(out(:,2)&out(:,3)==0,:)=[];
end

function r=CalculatedRatio(list_congru_active,conn_congru_active,list_congru_active_pair,congru_active_pair,FR,e)
    connidx_congru_active=list_congru_active(list_congru_active(:,2)>(FR-e)&list_congru_active(:,2)<(FR+e),:);
    pairidx_congru_active=list_congru_active_pair(list_congru_active_pair(:,2)>(FR-e)&list_congru_active_pair(:,2)<(FR+e));
    ratio=[];
    for j=1:114
        n1=findnumber([conn_congru_active(:,1);conn_congru_active(:,2)],connidx_congru_active(connidx_congru_active>=j*100000&connidx_congru_active(:,1)<(j+1)*100000));
        n2=findnumber([congru_active_pair(:,1);congru_active_pair(:,2)],pairidx_congru_active(pairidx_congru_active>=j*100000&pairidx_congru_active(:,1)<(j+1)*100000));
        ratio=[ratio;n1/2/n2];
    end
    ratio(ratio(:,1)==Inf|isnan(ratio(:,1))|ratio(:,1)==0,:)=[];
    r=ratio;
end    
    
function n=findnumber(b,s)
if isempty(s)
    n=0;
else
    b1=b(1:size(b,1)/2,1);
    b2=b(size(b,1)/2+1:end,1);
    temp=[];
    b1(:,2)=0;
    b2(:,2)=0;
    for i=1:size(s,1)
        b1(:,2)=double(b1(:,1)==s(i))+b1(:,2);
        b2(:,2)=double(b2(:,1)==s(i))+b2(:,2);
        %     temp=[temp,nnz(b1==s(i)& b2==s(i))];
    end
n=nnz(b1.*b2==1);
end
end

function conn_congru_active=FRchange(list_congru_active,conn_congru_active)
list_congru_active(:,4)=list_congru_active(:,2)-list_congru_active(:,3);
I_congru_active=list_congru_active(list_congru_active(:,4)>0,1);
D_congru_active=list_congru_active(list_congru_active(:,4)<0,1);

for i=1:size(conn_congru_active,1)
    conn_congru_active(i,4)=nnz(I_congru_active==conn_congru_active(i,1))*1+nnz(D_congru_active==conn_congru_active(i,1))*-1;
    conn_congru_active(i,5)=nnz(I_congru_active==conn_congru_active(i,2))*1+nnz(D_congru_active==conn_congru_active(i,2))*-1;
end
end