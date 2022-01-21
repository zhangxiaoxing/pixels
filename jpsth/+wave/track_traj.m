p3s=false;
mem_type='6s';

close all
data_str=wave.plot_cd('data_only',true,'mem_type',lower(mem_type));
fh=figure('Color','w','Position',[100,100,900,600]);
hold on;

proj1X3=data_str.proj1X3-mean([data_str.proj1X3;data_str.proj2X3]);
proj2X3=data_str.proj2X3-mean([data_str.proj1X3;data_str.proj2X3]);
proj1Y3=data_str.proj1Y3-mean([data_str.proj1Y3;data_str.proj2Y3]);
proj2Y3=data_str.proj2Y3-mean([data_str.proj1Y3;data_str.proj2Y3]);
proj1Z3=data_str.proj1Z3-mean([data_str.proj1Z3;data_str.proj2Z3]);
proj2Z3=data_str.proj2Z3-mean([data_str.proj1Z3;data_str.proj2Z3]);

proj1X6=data_str.proj1X6-mean([data_str.proj1X6;data_str.proj2X6]);
proj2X6=data_str.proj2X6-mean([data_str.proj1X6;data_str.proj2X6]);
proj1Y6=data_str.proj1Y6-mean([data_str.proj1Y6;data_str.proj2Y6]);
proj2Y6=data_str.proj2Y6-mean([data_str.proj1Y6;data_str.proj2Y6]);
proj1Z6=data_str.proj1Z6-mean([data_str.proj1Z6;data_str.proj2Z6]);
proj2Z6=data_str.proj2Z6-mean([data_str.proj1Z6;data_str.proj2Z6]);




timetag=1:0.03:44;
if p3s
    interpX2=pchip(1:44,proj2X3,1:0.03:44);
    interpY2=pchip(1:44,proj2Y3,1:0.03:44);
    interpZ2=pchip(1:44,proj2Z3,1:0.03:44);
    interpX1=pchip(1:44,proj1X3,1:0.03:44);
    interpY1=pchip(1:44,proj1Y3,1:0.03:44);
    interpZ1=pchip(1:44,proj1Z3,1:0.03:44);
else
    interpX2=pchip(1:44,proj2X6,1:0.03:44);
    interpY2=pchip(1:44,proj2Y6,1:0.03:44);
    interpZ2=pchip(1:44,proj2Z6,1:0.03:44);
    interpX1=pchip(1:44,proj1X6,1:0.03:44);
    interpY1=pchip(1:44,proj1Y6,1:0.03:44);
    interpZ1=pchip(1:44,proj1Z6,1:0.03:44);
end


if p3s
    v=VideoWriter(sprintf('track_traj_3s_s2_%s.mp4',mem_type),'MPEG-4');
else
    v=VideoWriter(sprintf('track_traj_6s_s2_%s.mp4',mem_type),'MPEG-4');
end
open(v);

for tt=1:numel(interpX1)
    cla
    grid on
    if p3s
        h13=plot3(interpX1(1:tt),interpY1(1:tt),interpZ1(1:tt),'-c','LineWidth',1.5);
        h23=plot3(interpX2(1:tt),interpY2(1:tt),interpZ2(1:tt),'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
        text(interpX1(tt),interpY1(tt),interpZ1(tt),'S1, 3s','Color','c','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top')
        text(interpX2(tt),interpY2(tt),interpZ2(tt),'S2, 3s','Color',[0.9290 0.6940 0.1250],'FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top')
    else
        h16=plot3(interpX1(1:tt),interpY1(1:tt),interpZ1(1:tt),'-k','LineWidth',1.5);
        h26=plot3(interpX2(1:tt),interpY2(1:tt),interpZ2(1:tt),'-r','LineWidth',1.5);
        text(interpX1(tt),interpY1(tt),interpZ1(tt),'S1, 6s','Color','k','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top')
        text(interpX2(tt),interpY2(tt),interpZ2(tt),'S2, 6s','Color','r','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top')
    end
    currt=timetag(tt);

    if currt<13
        tag='ITI';
        cc='k';
    elseif currt<17
        tag='Sample';
        cc='m';
    elseif currt<29
        tag='WM delay';
        cc='r';
    elseif currt<33
        tag='3s Test / 6s WM delay';
        if p3s
            cc='b';
        else
            cc='r';
        end
    elseif currt<41
        tag='6s WM delay';
        if p3s
            cc='k';
        else
            cc='r';
        end
    else
        tag='6s Test';
        if p3s
            cc='k';
        else
            cc='b';
        end
    end

    pdh=plot3(interpX1(tt),interpY1(tt),interpZ1(tt),'o','Color',cc,'MarkerSize',9);
    campos([interpX1(tt),interpY1(tt),interpZ1(tt)]*2)
    camtarget([0,0,0])
    
    if p3s
        switch(lower(mem_type))
            case 'both'
                xtickvec=-8:4:16;ytickvec=-8:4:16;ztickvec=-8:4:8;
            case '6s'
                xtickvec=-10:5:15;ytickvec=-10:5:10;ztickvec=-8:4:8;
            case '3s'
                xtickvec=-10:5:20;ytickvec=-8:4:12;ztickvec=-8:4:8;
        end
    else
        switch(lower(mem_type))
            case 'both'
                xtickvec=-5:5:15;ytickvec=-6:2:8;ztickvec=-8:4:8;
            case '6s'
                xtickvec=-2:1:2;ytickvec=-2:1:2;ztickvec=-8:4:8;
            case '3s'
                xtickvec=-4:5:16;ytickvec=-8:4:8;ztickvec=-8:4:8;
        end
    end
        set(gca(),'XTick',xtickvec,'YTick',ytickvec,'ZTick',ztickvec)
        xlim([min(xtickvec),max(xtickvec)])
        ylim([min(ytickvec),max(ytickvec)])
        zlim([min(ztickvec),max(ztickvec)])

    title(sprintf('%s selective, %s (%d msec)',mem_type,tag,int32(currt*250-4000)),'FontSize',20)
    pause(0.03)
    writeVideo(v,getframe(fh))
end
close(v);
