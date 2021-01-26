clearvars -Except FR_modulated
clc
close

%% extract licking information for ser files
f1=dir('D:\neupix\ser\*.ser');
lick_CR_matrix=cell(0);
lick_Hit_matrix=cell(0);
lick_matrix=cell(0);
for i=1:size(f1,1)
    try
    Data=ser2mat(fullfile('D:\neupix\ser',f1(i,1).name))';
    Sample=Data((Data(:,2)==9&Data(:,3)==1)|(Data(:,2)==10&Data(:,3)==100),:);
    Test=Data((Data(:,2)==9&Data(:,3)==3)|(Data(:,2)==10&Data(:,3)==2),:);
    Response=Data((Data(:,2)==4&Data(:,3)==1)|(Data(:,2)==5&Data(:,3)==1)|(Data(:,2)==6&Data(:,3)==1)|(Data(:,2)==7&Data(:,3)==1),2);
    Response(:,2)=0;
    Response(Response(:,1)==5|Response(:,1)==7,4)=1;
    a=40;
    while a<=length(Response)
        goodOff=nnz(Response(a-39:a,:)==5|Response(a-39:a,:)==7);
        if goodOff>=30 %.75 correct rate
            Response(a-39:a,2)=1;
        end
        a=a+1;
    end
    Response((Test-Sample)>5000,3)=1;
    
    Sample0=[];
    Response1=[];
    Sample0=Sample(Response(:,2)==1&Response(:,3)==1&Response(:,4)==1,:);   
    for t=1:size(Sample0)
        lick_matrix{end+1}=(Data((Data(:,1)>(Sample0(t,1)-3*1000))&(Data(:,1)<(Sample0(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample0(t,1));
    end
    Sample1=[];
    Response1=[];
    Sample1=Sample(Response(:,2)==1&Response(:,1)==5&Response(:,3)==1,:);
    Response1=Response(Response(:,2)==1&Response(:,1)==5&Response(:,3)==1,:);
%     lick_matrix=zeros(1,17*1000);
    for t=1:size(Sample1)
        lick_CR_matrix{end+1}=(Data((Data(:,1)>(Sample1(t,1)-3*1000))&(Data(:,1)<(Sample1(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample1(t,1));
    end
    lick=[];
    for t=1:size(Sample1)
        lick(t,:)=histcounts((Data((Data(:,1)>(Sample1(t,1)-3*1000))&(Data(:,1)<(Sample1(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample1(t,1))/1000,-3:0.2:14)/0.2;
    end    
    lick_CR(i,:)=mean(colfilt(lick,[1,5],'sliding',@(x)(sum(x)/5)),1);
    Sample2=[];
    Response2=[];
    Sample2=Sample(Response(:,2)==1&Response(:,1)==7&Response(:,3)==1,:);
    Response2=Response(Response(:,2)==1&Response(:,1)==7&Response(:,3)==1,:);
    for t=1:size(Sample2)
        lick_Hit_matrix{end+1}=(Data((Data(:,1)>(Sample2(t,1)-3*1000))&(Data(:,1)<(Sample2(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample2(t,1));
    end
    lick=[];
    for t=1:size(Sample2)
        lick(t,:)=histcounts((Data((Data(:,1)>(Sample2(t,1)-3*1000))&(Data(:,1)<(Sample2(t,1)+14*1000))&Data(:,2)==0&Data(:,3)==76,1)-Sample2(t,1))/1000,-3:0.2:14)/0.2;
    end
    lick_Hit(i,:)=mean(colfilt(lick,[1,5],'sliding',@(x)(sum(x)/5)),1);
    catch
        disp(f1(i,1).name)
        disp(i)
    end
end
save('lick_raster','lick_matrix')
save('lick6fromser','lick_CR','lick_Hit','lick_Hit_matrix','lick_CR_matrix')

%% 50 example trials of licking raster 
load('D:\code\lick_raster')
fh=figure('Color','w','Position',[100,100,500,500]);
a=0;b=0;
for i=1:50
%     if nnz(lick_matrix{i}/10>800 & lick_matrix{i}/10<900)>0
%         a=a+1;
%     else
%         b=b+1;
%     end
    if nnz(lick_matrix{i}/10>800 & lick_matrix{i}/10<900)>0
        a =fill([-300:1400;-300:1400]',[50+1-i-0.5,50+1-i+0.5], 'r','edgecolor',[0.98 0.9 0.9]);
        hold on
        plot([lick_matrix{i}/10 lick_matrix{i}/10], [50+1-i-0.5 50+1-i+0.5],'k');
        hold on
    else
        a =fill([-300:1400;-300:1400]',[50+1-i-0.5,50+1-i+0.5], 'r','edgecolor',[0.8 0.9 0.98]);
        try
            plot([lick_matrix{i}/10 lick_matrix{i}/10], [50+1-i-0.5 50+1-i+0.5],'k');
            hold on
        catch
        end
    end
end
plot([0 0],[0 100],'k--')
plot([100 100],[0 100],'k--')
plot([700 700],[0 100],'k--')
plot([800 800],[0 100],'k--')
plot([900 900],[0 100],'k--')
xlim([-300 1400])
ylim([0 50])
set(gca,'XTick',-300:100:1400,'XTickLabel',{' ',' ',' ','0',' ',' ',' '' ',' ',' ','6',' ',' ',' ',' '' ',' ','12',' ',' '})
xlabel('Time (s)')
ylabel('Trial ID #')
exportgraphics(fh,'lick_raster.pdf','ContentType','vector');

%% licking rate of all sessions
load('lick6fromser.mat')
fh=figure('Color','w','Position',[100,100,500,500]);
hold on
time=size(lick_Hit,2);
Time = [1:time, fliplr(1:time)];
Highervalue = smooth(mean(lick_Hit,1))' + std(lick_Hit,0,1)/sqrt(size(lick_Hit,1));
Lowervalue = smooth(mean(lick_Hit,1))' - std(lick_Hit,0,1)/sqrt(size(lick_Hit,1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, 'r', 'edgecolor','none');
alpha(a,0.2);
plot(1:time,smooth(mean(lick_Hit,1))','r-')
hold on
time=size(lick_CR,2);
Time = [1:time, fliplr(1:time)];
Highervalue = smooth(mean(lick_CR,1))' + std(lick_CR,0,1)/sqrt(size(lick_CR,1));
Lowervalue = smooth(mean(lick_CR,1))' - std(lick_CR,0,1)/sqrt(size(lick_CR,1));
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, 'b', 'edgecolor','none');
alpha(a,0.2);
plot(1:time,smooth(mean(lick_CR,1))','b-')

plot([15 15],[0 8],'k--')
plot([20 20],[0 8],'k--')
plot([35 35],[0 8],'k--')
plot([40 40],[0 8],'k--')
plot([45 45],[0 8],'k--')
plot([50 50],[0 8],'k--')
plot([55 55],[0 8],'k--')
plot([60 60],[0 8],'k--')
set(gca,'XTick',1:5:85,'XTickLabel',{' ',' ',' ','0',' ',' ',' ',' ',' ','6',' ',' ',' ','','','12'})
xlabel('Time(s)','FontSize',10);
ylabel('Lick(Hz)','FontSize',10);
box off
exportgraphics(fh,'lick_rate.pdf','ContentType','vector');

