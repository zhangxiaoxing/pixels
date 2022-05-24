clc
clear
CD=cd;
foldername=dir('*');
% foldername='D:\neupixel\DataSum\191021-DPA-Combined_Learning8_31_g0';
Perf=[];
path_list=cell(0);
foldername(startsWith({foldername.name},'.') | ~[foldername.isdir] | startsWith({foldername.name},'counts.npz') | startsWith({foldername.name},'sync.ffs_db'))=[];
for folder_ID=1:size(foldername,1)
    
    folder1name=dir([CD '\' foldername(folder_ID).name '\*cleaned']);
    folder1name(startsWith({folder1name.name},'.') | ~[folder1name.isdir])=[];
    
    CD2=[folder1name(1).folder '\' folder1name(1).name '\'];
    %%    
    trials = double(hdf5read([CD2 'events.hdf5'],'/trials')) ;
    trials(9,:)=trials(5,:)-trials(6,:);
    trials(9,trials(9,:)==0)=-1;
    trials(9,trials(9,:)==4|trials(9,:)==-4)=1;
    trials(10,:)=0;
    trials(10,trials(7,:)==trials(9,:))=1;
    trials(11,:)=0;
    a=40;
    while a<=length(trials)
        goodOff=nnz(xor(trials(5,a-39:a)==trials(6,a-39:a),trials(7,a-39:a)>0));
        if goodOff>=30 %.75 correct rate
            trials(11,a-39:a)=1;
        end
        a=a+1;
    end
   
    Tiralnum_3=length(find(trials(8,:)==3&trials(11,:)==1));
    Tiralnum_6=length(find(trials(8,:)==6&trials(11,:)==1));
    if Tiralnum_3==0 || Tiralnum_6==0
    else
        Perf_3(:,1)=length(find(trials(8,:)==3&trials(10,:)==1&trials(11,:)==1))/Tiralnum_3;
        Perf_6(:,1)=length(find(trials(8,:)==6&trials(10,:)==1&trials(11,:)==1))/Tiralnum_6;
        Perf_3(:,2)=length(find(trials(8,:)==3&trials(9,:)==1&trials(10,:)==1&trials(11,:)==1))/length(find(trials(8,:)==3&trials(9,:)==1&trials(11,:)==1));
        Perf_6(:,2)=length(find(trials(8,:)==6&trials(9,:)==1&trials(10,:)==1&trials(11,:)==1))/length(find(trials(8,:)==6&trials(9,:)==1&trials(11,:)==1));
        Perf_3(:,3)=length(find(trials(8,:)==3&trials(9,:)==-1&trials(10,:)==1&trials(11,:)==1))/length(find(trials(8,:)==3&trials(9,:)==-1&trials(11,:)==1));
        Perf_6(:,3)=length(find(trials(8,:)==6&trials(9,:)==-1&trials(10,:)==1&trials(11,:)==1))/length(find(trials(8,:)==6&trials(9,:)==-1&trials(11,:)==1));
        Perf=[Perf;[Perf_3,Perf_6]];
        path_list{end+1,1}=CD2;
%         Perf{folder_ID,1}=foldername(folder_ID).name;
%         Perf{folder_ID,2}=[Perf_3,Perf_6];
    end
    clearvars -Except Perf foldername folder_ID CD path_list
end
% bar(mean(Perf,1))
% hold on
mouseID=unique(Perf(:,3));
for i=1:size(mouseID,1)
    Perf_3(i,:)=mean(Perf(Perf(:,7)==mouseID(i,1),1),1);
    Perf_6(i,:)=mean(Perf(Perf(:,7)==mouseID(i,1),2),1);
%     plot([Perf_3(i,:), Perf_6(i,:)],'color','k','marker','.','MarkerSize',2,'linewidth',1)
%     hold on
end
% p=signrank(Perf(:,1),Perf(:,2))
% errorbar(1:2,mean(Perf,1),std(Perf)/sqrt(size(Perf,1)),'LineStyle','none','color','k','marker','o','MarkerSize',5,'linewidth',2)

xlim([0,3])
ylim([0.5 1])
set(gca,'XTick',1:2,'XTickLabel',{'3s','6s'})
ylabel('Performance (%)')
title('Performance n=40')

saveas(gcf,'D:\code','fig')
