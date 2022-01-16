p3s=true;
close all
load viewport_data
fh=figure('Color','w','Position',[100,100,900,600]);
hold on;

proj1X3=proj1X3-mean([proj1X3;proj2X3]);
proj2X3=proj2X3-mean([proj1X3;proj2X3]);
proj1Y3=proj1Y3-mean([proj1Y3;proj2Y3]);
proj2Y3=proj2Y3-mean([proj1Y3;proj2Y3]);
proj1Z3=proj1Z3-mean([proj1Z3;proj2Z3]);
proj2Z3=proj2Z3-mean([proj1Z3;proj2Z3]);


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
    v=VideoWriter('track_traj_3s_s2.mp4','MPEG-4');
else
    v=VideoWriter('track_traj_6s_s2.mp4','MPEG-4');
end
open(v);

for tt=1:numel(interpZ)
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

    pdh=plot3(interpX2(tt),interpY2(tt),interpZ2(tt),'o','Color',cc,'MarkerSize',9);
    campos([interpX2(tt),interpY2(tt),interpZ2(tt)]*2)
    camtarget([0,0,0])
    currt=timetag(tt);
    if p3s
        set(gca(),'XTick',-8:4:16,'YTick',-8:4:16,'ZTick',-8:4:8)
        xlim([-8,16])
        ylim([-8,16])
        zlim([-8,8])
    else
        set(gca(),'XTick',-5:5:15,'YTick',-6:2:8,'ZTick',-6:2:6)
        xlim([-5,15])
        ylim([-6,8])
        zlim([-6,6])
    end
    title(sprintf('Both selective, %s (%d msec)',tag,int32(currt*250-4000)),'FontSize',20)
    pause(0.03)
    writeVideo(v,getframe(fh))
end
close(v);
