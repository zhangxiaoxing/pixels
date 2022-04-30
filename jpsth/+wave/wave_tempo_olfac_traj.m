% FR->D3S1 D6S1 D3S2 D6S2
function wave_tempo_olfac_traj(opt)
arguments
    opt.pca (1,1) logical = false
    opt.gen_movie (1,1) logical = false
end


normmat=normalize(pcamat.','range');
if opt.pca
    %% PCA
    [coef,score,latent]=pca(normmat);
    fh=figure('Color','w');
    hold on
    d3s1h=plot3(score(1:56,1),score(1:56,2),score(1:56,3),'--r');
    d3s2h=plot3(score((1:56)+56,1),score((1:56)+56,2),score((1:56)+56,3),'--b');
    d6s1h=plot3(score((1:56)+112,1),score((1:56)+112,2),score((1:56)+112,3),'-r');
    d6s2h=plot3(score((1:56)+168,1),score((1:56)+168,2),score((1:56)+168,3),'-b');
    for tt=2:2:14
        text(score(tt*4,1),score(tt*4,2),score(tt*4,3)-4,num2str(tt-4),'FontSize',12);
    end
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    legend([d6s1h,d6s2h,d3s1h,d3s2h],{'S1 6s','S2 6s','S1 3s','S2 3s'},Location='northoutside',Orientation='horizontal')
    view(20,10)
    if opt.gen_movie
        v=VideoWriter('tempo_olfac_pca.mp4','MPEG-4');
        open(v)
        for az=20:0.5:380
            view(az,10)
            writeVideo(v,getframe(fh))
        end
        close(v)
    end
end

%% coding direction
% odor cd
ocd=mean(normmat([17:28,(17:40)+112],:))-mean(normmat([(17:28)+56,(17:40)+168],:));
ocd=ocd/norm(ocd);
% duration cd
dcd=mean(normmat([1:12,(1:12)+56],:))-mean(normmat([(1:12)+112,(1:12)+168],:));
dcd=dcd/norm(dcd);
if opt.gen_movie
    v=VideoWriter('tempo_olfac_cd_develop.mp4','MPEG-4');
    open(v)
end
fh=figure('Color','w','Position',[32,32,960,400]);
for ep=1:56
    for subidx=1:2
        subplot(1,2,subidx)
        cla
        hold on
        grid on
        d3s1h=plot3(normmat(1:ep,:)*(ocd.'),normmat(1:ep,:)*(dcd.'),1:ep,'--r');
        d3s2h=plot3(normmat((1:ep)+56,:)*(ocd.'),normmat((1:ep)+56,:)*(dcd.'),1:ep,'--b');
        d6s1h=plot3(normmat((1:ep)+112,:)*(ocd.'),normmat((1:ep)+112,:)*(dcd.'),1:ep,'-r');
        d6s2h=plot3(normmat((1:ep)+168,:)*(ocd.'),normmat((1:ep)+168,:)*(dcd.'),1:ep,'-b');

        %     text(normmat(20,:)*(ocd.'),normmat(20,:)*(dcd.'),20,'S1,3s','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle');
        %     text(normmat(24+56,:)*(ocd.'),normmat(24+56,:)*(dcd.'),24,'S2,3s','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle');
        %     text(normmat(28+112,:)*(ocd.'),normmat(28+112,:)*(dcd.'),28,'S1,6s','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle');
        %     text(normmat(32+168,:)*(ocd.'),normmat(32+168,:)*(dcd.'),32,'S2,6s','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle');

        set(gca,'ZTick',8:8:56,'ZTickLabel',-2:2:10)
        zlabel('Time (s)','FontSize',16)
        if subidx==1
            view([185,5])
            xlabel('Olfactory CD','FontSize',16)
            ylabel('')
        else
            xlabel('')
            ylabel('Duration prediction CD','FontSize',16)
            view([275,5])
            if ep<12
                sgtitle('ITI','FontSize',20,'FontWeight','bold')
            elseif ep<16
                sgtitle('Sample','FontSize',20,'FontWeight','bold')
            elseif ep<28
                sgtitle('Delay','FontSize',20,'FontWeight','bold')
            elseif ep<32
                sgtitle('3s: Test, 6s: Delay','FontSize',20,'FontWeight','bold')
            elseif ep<40
                sgtitle('6s: Delay','FontSize',20,'FontWeight','bold')
            elseif ep<44
                sgtitle('6s: Test','FontSize',20,'FontWeight','bold')
            else
                sgtitle('')
            end
        end
    end
    pause(0.2)
    if opt.gen_movie
        for frame=1:12
            writeVideo(v,getframe(fh))
        end
    end
end
% legend([d6s1h,d6s2h,d3s1h,d3s2h],{'S1 6s','S2 6s','S1 3s','S2 3s'},Location='northoutside',Orientation='horizontal')
if opt.gen_movie

    close(v)

    v=VideoWriter('tempo_olfac_cd.mp4','MPEG-4');
    open(v)
    for az=-30:0.5:330
        view(az,15)
        writeVideo(v,getframe(fh))
    end
    close(v)
end
end




