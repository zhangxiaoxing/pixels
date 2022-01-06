function plot_cd(opt)
arguments
    opt.plot_rpt (1,1) double {mustBePositive,mustBeInteger} =2
    opt.stats_rpt (1,1) double {mustBePositive,mustBeInteger} =5
    opt.plot_traj (1,1) logical = true
    opt.alt_proj (1,1) logical = false
end
% persistent com_str onepath_ delay_ selidx_ decision_ rnd_half_ curve_


% meta3=ephys.util.load_meta('type','neupix','delay',3);
% meta6=ephys.util.load_meta('type','neupix','delay',6);
[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir('type','raw');

[s1mat3,s2mat3,s1mat6,s2mat6,s1mat31,s2mat31,s1mat61,s2mat61,waveids]=deal([]);
[s1shufmat,s2shufmat]=deal(cell(1,2*opt.stats_rpt));

for ii=reshape(cell2mat(sessmap.keys()),1,[])
    disp(ii)
    fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
    fr=h5read(fpath,'/FR_All');
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    dur_resp=behav.tag_block(trials,'wt',false);

    wavesess=ones(size(suid))*double(ii);
    waveidsess=ephys.get_wave_id(wavesess,suid,'early',false,'ctx',true);
%     s1sel3=waveids==1;
%     s1sel6=waveids==3;
%     s1selb=waveids==5;
%     
%     s2sel3=waveids==2;
%     s2sel6=waveids==4;
%     s2selb=waveids==6;


    s1d3=trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0 & dur_resp(:,3)>1;
    s2d3=trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0 & dur_resp(:,3)>1;
    s1d6=trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0 & dur_resp(:,3)>1;
    s2d6=trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0 & dur_resp(:,3)>1;

    s1d31=trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & dur_resp(:,3)==1;
    s2d31=trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & dur_resp(:,3)==1;
    s1d61=trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & dur_resp(:,3)==1;
    s2d61=trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & dur_resp(:,3)==1;

    if false
        for rpt=1:opt.stats_rpt
            s1pool=find(s1d3);
            s2pool=find(s2d3);
            s1a=randsample(s1pool,floor(numel(s1pool)./2));
            s1b=s1pool(~ismember(s1pool,s1a));
            s2a=randsample(s2pool,floor(numel(s2pool)./2));
            s2b=s2pool(~ismember(s2pool,s2a));
            s1shufmat{rpt*2-1}=[s1shufmat{rpt*2-1};shiftdim(mean(fr(s1a,msel,1:44)),1)];
            s1shufmat{rpt*2}=[s1shufmat{rpt*2};shiftdim(mean(fr(s1b,msel,1:44)),1)];
            s2shufmat{rpt*2-1}=[s2shufmat{rpt*2-1};shiftdim(mean(fr(s2a,msel,1:44)),1)];
            s2shufmat{rpt*2}=[s2shufmat{rpt*2};shiftdim(mean(fr(s2b,msel,1:44)),1)];
        end
    end
    wavesel=ismember(waveidsess,1:6);
    if ~any(wavesel)
        continue;
    end

    s1mat3=[s1mat3;shiftdim(mean(fr(s1d3,wavesel,1:44),1),1)];
    s2mat3=[s2mat3;shiftdim(mean(fr(s2d3,wavesel,1:44),1),1)];
    s1mat6=[s1mat6;shiftdim(mean(fr(s1d6,wavesel,1:44),1),1)];
    s2mat6=[s2mat6;shiftdim(mean(fr(s2d6,wavesel,1:44),1),1)];

    s1mat31=[s1mat31;shiftdim(mean(fr(s1d31,wavesel,1:44),1),1)];
    s2mat31=[s2mat31;shiftdim(mean(fr(s2d31,wavesel,1:44),1),1)];
    s1mat61=[s1mat61;shiftdim(mean(fr(s1d61,wavesel,1:44),1),1)];
    s2mat61=[s2mat61;shiftdim(mean(fr(s2d61,wavesel,1:44),1),1)];


    if ~(isequal(size(s1mat3),size(s1mat31)) ...
            && isequal(size(s2mat3),size(s2mat31))...
            && isequal(size(s1mat6),size(s1mat61))...
            && isequal(size(s2mat6),size(s2mat61)))
        keyboard()
    end

    waveids=[waveids;waveidsess(wavesel)];

end

normfr=normalize([s1mat3,s2mat3,s1mat6,s2mat6,s1mat31,s2mat31,s1mat61,s2mat61].','range').';
[s1mat3,s2mat3,s1mat6,s2mat6]=deal(normfr(:,1:44),normfr(:,(1:44)+44),normfr(:,(1:44)+88),normfr(:,(1:44)+132));
[s1mat31,s2mat31,s1mat61,s2mat61]=deal(normfr(:,(1:44)+44*4),normfr(:,(1:44)+44*5),normfr(:,(1:44)+44*6),normfr(:,(1:44)+44*7));
if opt.alt_proj
    %3s vs 6s axis
else
    cdMat=[s2mat3,s2mat6]-[s1mat3,s1mat6];
    cdDelay=mean(cdMat(:,[17:28,(17:40)+44]),2);
    cdDelay=cdDelay/norm(cdDelay);
    [cdDelay3,cdDelay6,cdDelayB]=deal(cdDelay);
    cdDelay3(~ismember(waveids,1:2))=0;
    cdDelay6(~ismember(waveids,3:4))=0;
    cdDelayB(~ismember(waveids,5:6))=0;
    
%     cdDelayM=mean(cdMat(:,21:24),2);
%     cdDelayM=cdDelayM/norm(cdDelayM);
%     
%     cdDelayL=mean(cdMat(:,25:28),2);
%     cdDelayL=cdDelayL/norm(cdDelayL);
end
if opt.plot_traj
    if false
        for rpt=1:(opt.plot_rpt*2)
            [proj1E,proj1M,proj1L]=deal(s1shufmat{rpt}.'*cdDelay,s1shufmat{rpt}.'*cdDelayM,s1shufmat{rpt}.'*cdDelayL);
            [proj2E,proj2M,proj2L]=deal(s2shufmat{rpt}.'*cdDelay,s2shufmat{rpt}.'*cdDelayM,s2shufmat{rpt}.'*cdDelayL);
            %     plot3(proj1E(1:5),proj1M(1:5),proj1L(1:5),':r.','LineWidth',0.5)
            %     plot3(proj2E(1:5),proj2M(1:5),proj2L(1:5),':b.','LineWidth',0.5)
            %     plot3(proj1E(5:9),proj1M(5:9),proj1L(5:9),'--r.','LineWidth',0.5)
            %     plot3(proj2E(5:9),proj2M(5:9),proj2L(5:9),'--b.','LineWidth',0.5)
            plot3(smooth(proj1E,3),smooth(proj1M,3),smooth(proj1L,3),'-m','LineWidth',0.5)
            plot3(smooth(proj2E,3),smooth(proj2M,3),smooth(proj2L,3),'-c','LineWidth',0.5)
            %     plot3(proj1E(9),proj1M(9),proj1L(9),'ko','MarkerFaceColor','k')
            %     plot3(proj2E(9),proj2M(9),proj2L(9),'ko','MarkerFaceColor','k')
        end
    end
    % plot3(proj1E(1:5),proj1M(1:5),proj1L(1:5),':r.','LineWidth',1)
    % plot3(proj2E(1:5),proj2M(1:5),proj2L(1:5),':b.','LineWidth',1)

    [proj1E3,proj1M3,proj1L3]=deal(smooth(s1mat3.'*cdDelay3),smooth(s1mat3.'*cdDelay6),smooth(s1mat3.'*cdDelayB));
    [proj2E3,proj2M3,proj2L3]=deal(smooth(s2mat3.'*cdDelay3),smooth(s2mat3.'*cdDelay6),smooth(s2mat3.'*cdDelayB));

    [proj1E6,proj1M6,proj1L6]=deal(smooth(s1mat6.'*cdDelay3),smooth(s1mat6.'*cdDelay6),smooth(s1mat6.'*cdDelayB));
    [proj2E6,proj2M6,proj2L6]=deal(smooth(s2mat6.'*cdDelay3),smooth(s2mat6.'*cdDelay6),smooth(s2mat6.'*cdDelayB));

    [proj1E31,proj1M31,proj1L31]=deal(smooth(s1mat31.'*cdDelay3),smooth(s1mat31.'*cdDelay6),smooth(s1mat31.'*cdDelayB));
    [proj2E31,proj2M31,proj2L31]=deal(smooth(s2mat31.'*cdDelay3),smooth(s2mat31.'*cdDelay6),smooth(s2mat31.'*cdDelayB));
    [proj1E61,proj1M61,proj1L61]=deal(smooth(s1mat61.'*cdDelay3),smooth(s1mat61.'*cdDelay6),smooth(s1mat61.'*cdDelayB));
    [proj2E61,proj2M61,proj2L61]=deal(smooth(s2mat61.'*cdDelay3),smooth(s2mat61.'*cdDelay6),smooth(s2mat61.'*cdDelayB));




    v=VideoWriter('cd_devp_movie.mp4','MPEG-4');
    open(v);
    fh=figure('Color','w','Position',[100,100,900,600]);
    hold on;
    for tt=44%1:44
        cla
        view([-5,10])
        h13=plot3(proj1E3(1:tt),proj1M3(1:tt),proj1L3(1:tt),'-c','LineWidth',1.5);
        h23=plot3(proj2E3(1:tt),proj2M3(1:tt),proj2L3(1:tt),'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
        h16=plot3(proj1E6(1:tt),proj1M6(1:tt),proj1L6(1:tt),'-k','LineWidth',1.5);
        h26=plot3(proj2E6(1:tt),proj2M6(1:tt),proj2L6(1:tt),'-r','LineWidth',1.5);
        legend([h26,h23,h16,h13],{'S1-6s trial','S1-3s trial','S2-6s trial','S2-3s trial'},'FontSize',16,'Location','eastoutside','Orientation','vertical','AutoUpdate','off')

        th1=text(proj1E3(tt),proj1M3(tt),proj1L3(tt),'S2-3s','Color','c','FontSize',16);
        th2=text(proj2E3(tt),proj2M3(tt),proj2L3(tt),'S1-3s','Color',[0.9290 0.6940 0.1250],'FontSize',16);
        th3=text(proj1E6(tt),proj1M6(tt),proj1L6(tt),'S2-6s','Color','k','FontSize',16);
        th4=text(proj2E6(tt),proj2M6(tt),proj2L6(tt),'S1-6s','Color','r','FontSize',16);
        pause(0.5)
        if tt<13
            tag='ITI';
        elseif tt<17
            tag='Sample';
        elseif tt<29
            tag='WM delay';
        elseif tt<33
            tag='3s Test / 6s WM delay';
        elseif tt<41
            tag='6s WM delay';
        else
            tag='6s Test';
        end
        title(sprintf('%s (%d msec)',tag,tt*250-4000),'FontSize',20)
        xlabel('3s selective')
        ylabel('6s selective')
        zlabel('Both selective')
        for frame=1:15
%             writeVideo(v,getframe(fh))
        end
    end
    close(v);
    keyboard

    delete(th1)
    delete(th2)
    delete(th3)
    delete(th4)

    th1=text(proj1E3(20),proj1M3(20),proj1L3(20),'S2-3s','Color','c','FontSize',16);
    th2=text(proj2E3(24),proj2M3(24),proj2L3(24),'S1-3s','Color',[0.9290 0.6940 0.1250],'FontSize',16);
    th3=text(proj1E6(28),proj1M6(28),proj1L6(28),'S2-6s','Color','k','FontSize',16);
    th4=text(proj2E6(32),proj2M6(32),proj2L6(32),'S1-6s','Color','r','FontSize',16);

    v=VideoWriter('cd_rotate_movie.mp4','MPEG-4');
    open(v);
    for dg=-5:10:355
        view([dg,10])
        for fm=1:10
            writeVideo(v,getframe(fh))
        end
    end
    close(v)


    h131=plot3(proj1E31(1:tt),proj1M31(1:tt),proj1L31(1:tt),'--c','LineWidth',1.5);
    h231=plot3(proj2E31(1:tt),proj2M31(1:tt),proj2L31(1:tt),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
    h161=plot3(proj1E61(1:tt),proj1M61(1:tt),proj1L61(1:tt),'--k','LineWidth',1.5);
    h261=plot3(proj2E61(1:tt),proj2M61(1:tt),proj2L61(1:tt),'--r','LineWidth',1.5);

    v=VideoWriter('cd_rotate_movie_blockN1.mp4','MPEG-4');
    open(v);
    for dg=-5:10:355
        view([dg,10])
        for fm=1:10
            writeVideo(v,getframe(fh))
        end
    end
    close(v)




    distmat=[proj1E6-proj1E3,proj1M6-proj1M3,proj1L6-proj1L3];
    distmat1=[proj1E61-proj1E31,proj1M61-proj1M31,proj1L61-proj1L31];
    dist_3_6=arrayfun(@(x) norm(distmat(x,:)),1:44);
    dist_31_61=arrayfun(@(x) norm(distmat1(x,:)),1:44);

    distmat2=[proj2E6-proj2E3,proj2M6-proj2M3,proj2L6-proj2L3];
    distmat21=[proj2E61-proj2E31,proj2M61-proj2M31,proj2L61-proj2L31];
    dist_3_6_2=arrayfun(@(x) norm(distmat2(x,:)),1:44);
    dist_31_61_2=arrayfun(@(x) norm(distmat21(x,:)),1:44);

    f36=figure('Color','w');
    hold on
    ph1=plot(dist_3_6,'r-');
    ph11=plot(dist_31_61,'r--');
    ph2=plot(dist_3_6_2,'b-');
    ph21=plot(dist_31_61_2,'b--');
    ylabel('CD trajectory distance (AU)')
    set(gca,'XTick',16.5:20:45,'XTickLabel',[0,5])
    xlabel('Time (s)')
    arrayfun(@(x) xline(x,':k'),[12.5,16.5,28.5,32.5,40.5,44.5])
    legend([ph1,ph2,ph11,ph21],{'S1 3s v.s. S1 6s','S2 3s v.s. S2 6s','S1 3/6s, 1st trial in block','S2 3/6s, 1st trial in block'},'Location','eastoutside')
    exportgraphics(f36,'cd_dist_3s_6s.pdf','ContentType','vector')


    distmat=[proj1E6-proj2E6,proj1M6-proj2M6,proj1L6-proj2L6];
    distmat1=[proj1E61-proj2E61,proj1M61-proj2M61,proj1L61-proj2L61];
    dist_3_6=arrayfun(@(x) norm(distmat(x,:)),1:44);
    dist_31_61=arrayfun(@(x) norm(distmat1(x,:)),1:44);

    distmat2=[proj2E6-proj2E3,proj2M6-proj2M3,proj2L6-proj2L3];
    distmat21=[proj2E61-proj2E31,proj2M61-proj2M31,proj2L61-proj2L31];
    dist_3_6_2=arrayfun(@(x) norm(distmat2(x,:)),1:44);
    dist_31_61_2=arrayfun(@(x) norm(distmat21(x,:)),1:44);

    f36=figure('Color','w');
    hold on
    ph1=plot(dist_3_6,'r-');
    ph11=plot(dist_31_61,'r--');
    ph2=plot(dist_3_6_2,'b-');
    ph21=plot(dist_31_61_2,'b--');
    ylabel('CD trajectory distance (AU)')
    set(gca,'XTick',16.5:20:45,'XTickLabel',[0,5])
    xlabel('Time (s)')
    arrayfun(@(x) xline(x,':k'),[12.5,16.5,28.5,32.5,40.5,44.5])
    legend([ph1,ph2,ph11,ph21],{'S1 3s v.s. S1 6s','S2 3s v.s. S2 6s','S1 3/6s, 1st trial in block','S2 3/6s, 1st trial in block'},'Location','eastoutside')
    exportgraphics(f36,'cd_dist_3s_6s.pdf','ContentType','vector')


%     plot3(proj1E(1:9),proj1M(1:9),proj1L(1:9),':r.','LineWidth',1)
%     plot3(proj2E(1:9),proj2M(1:9),proj2L(1:9),':b.','LineWidth',1)
%     plot3(proj1E(9:end),proj1M(9:end),proj1L(9:end),'-r.','LineWidth',1)
%     plot3(proj2E(9:end),proj2M(9:end),proj2L(9:end),'-b.','LineWidth',1)
% 
%     plot3(proj1E(9),proj1M(9),proj1L(9),'ko','MarkerFaceColor','k')
%     plot3(proj2E(9),proj2M(9),proj2L(9),'ko','MarkerFaceColor','k')

%     xlabel('Early CD (A.U.)')
%     ylabel('Mid CD (A.U.)')
%     zlabel('Late CD (A.U.)')
%     grid on
%     keyboard()
%     exportgraphics(fh,'CD_project_traj_3d.pdf','ContentType','vector');
end
estat=[];
mstat=[];
lstat=[];

for rpt=1:opt.stats_rpt
    [proj1E,proj1M,proj1L]=deal(s1shufmat{rpt}.'*cdDelay,s1shufmat{rpt}.'*cdDelayM,s1shufmat{rpt}.'*cdDelayL);
    [proj2E,proj2M,proj2L]=deal(s2shufmat{rpt}.'*cdDelay,s2shufmat{rpt}.'*cdDelayM,s2shufmat{rpt}.'*cdDelayL);
    estat=[estat,proj1E-proj2E];
    mstat=[mstat,proj1M-proj2M];
    lstat=[lstat,proj1L-proj2L];
end
[emm,mmm,lmm]=deal(mean(estat,2),mean(mstat,2),mean(lstat,2));
eci=prctile(estat.',[2.5,97.5]);
mci=prctile(mstat.',[2.5,97.5]);
lci=prctile(lstat.',[2.5,97.5]);

fh=figure('Color','w','Position',[100,100,190,150]);
hold on
ph=plot([emm,mmm,lmm]);
[ph(1).Color,ph(2).Color,ph(3).Color]=deal('r','b','k');
fill([1:32,32:-1:1],[eci(1,:),fliplr(eci(2,:))],'r','EdgeColor','none','FaceAlpha',0.1)
fill([1:32,32:-1:1],[mci(1,:),fliplr(mci(2,:))],'b','EdgeColor','none','FaceAlpha',0.1)
fill([1:32,32:-1:1],[lci(1,:),fliplr(lci(2,:))],'k','EdgeColor','none','FaceAlpha',0.1)
set(gca(),'XTick',8.5:20:32,'XTickLabel',0:5:5)
xline(4.5,'--k');xline(8.5,'--k')
xlabel('Time (s)')
ylabel('CD distance (A.U.)')
exportgraphics(fh,'CD_project_dist.pdf','ContentType','vector');
keyboard();


end
