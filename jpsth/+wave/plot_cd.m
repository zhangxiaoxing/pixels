function data_str=plot_cd(opt)
arguments
    opt.plot_rpt (1,1) double {mustBePositive,mustBeInteger} =2
    opt.stats_rpt (1,1) double {mustBePositive,mustBeInteger} =5
    opt.plot_traj (1,1) logical = true
    opt.alt_proj (1,1) logical = false
    opt.cd_pc_proj (1,1) logical = false
    opt.cd_tcom_proj (1,1) logical = true
    opt.plot_1st_trial (1,1) logical = false
    opt.gen_movie (1,1) logical = false
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'mem','3s','6s','both'})} = 'both'
    opt.data_only (1,1) logical = false
    
end
% persistent com_str onepath_ delay_ selidx_ decision_ rnd_half_ curve_
if nargout>0
    data_str=struct();
end

[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir('type','raw');

[s1_FR_3,s2_FR_3,s1_FR_6,s2_FR_6,s1_FR_3_t1,s2_FR_3_t1,s1_FR_6_t1,s2_FR_6_t1,waveids,comalls1,comalls2]=deal([]);
[s1shufmat,s2shufmat]=deal(cell(1,2*opt.stats_rpt));

if opt.cd_tcom_proj
    commap=wave.get_com_map();
end

for ii=reshape(cell2mat(sessmap.keys()),1,[])
    
    if opt.cd_tcom_proj && ~isfield(commap,sprintf('s%d',ii))
        continue
    end
    disp(ii)
    s1commap=commap.(sprintf('s%d',ii)).s1;
    s2commap=commap.(sprintf('s%d',ii)).s2;

    fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
    fr=h5read(fpath,'/FR_All');
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    dur_resp=behav.tag_block(trials,'wt',false);

    wavesess=ones(size(suid))*double(ii);
    waveidsess=ephys.get_wave_id(wavesess,suid,'early',false,'ctx',true);
    if opt.plot_1st_trial
    s1d3=trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0 & dur_resp(:,3)>1;
    s2d3=trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0 & dur_resp(:,3)>1;
    s1d6=trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0 & dur_resp(:,3)>1;
    s2d6=trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0 & dur_resp(:,3)>1;

    s1d31=trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & dur_resp(:,3)==1;
    s2d31=trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & dur_resp(:,3)==1;
    s1d61=trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & dur_resp(:,3)==1;
    s2d61=trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & dur_resp(:,3)==1;
    else
        s1d3=trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0;
        s2d3=trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0;
        s1d6=trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0;
        s2d6=trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0;
    end
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
    switch opt.mem_type
        case 'mem'
            wavesel=ismember(waveidsess,1:6);
        case '3s'
            wavesel=ismember(waveidsess,1:2);
        case '6s'
            wavesel=ismember(waveidsess,3:4);
        case 'both'
            wavesel=ismember(waveidsess,5:6);
    end

    if ~any(wavesel)
        continue;
    end

    s1_FR_3=[s1_FR_3;shiftdim(mean(fr(s1d3,wavesel,1:44),1),1)];
    s2_FR_3=[s2_FR_3;shiftdim(mean(fr(s2d3,wavesel,1:44),1),1)];
    s1_FR_6=[s1_FR_6;shiftdim(mean(fr(s1d6,wavesel,1:44),1),1)];
    s2_FR_6=[s2_FR_6;shiftdim(mean(fr(s2d6,wavesel,1:44),1),1)];
    if opt.cd_tcom_proj
        currid=suid(wavesel);
        coms1=zeros(nnz(wavesel),1);
        coms2=zeros(nnz(wavesel),1);
        coms1(s1commap.isKey(num2cell(currid)))=cell2mat(s1commap.values(num2cell(currid(isKey(s1commap,num2cell(currid))))));
        coms2(s2commap.isKey(num2cell(currid)))=cell2mat(s2commap.values(num2cell(currid(isKey(s2commap,num2cell(currid))))));
        comalls1=[comalls1;coms1];
        comalls2=[comalls2;coms2];
    end

    if opt.plot_1st_trial
        s1_FR_3_t1=[s1_FR_3_t1;shiftdim(mean(fr(s1d31,wavesel,1:44),1),1)];
        s2_FR_3_t1=[s2_FR_3_t1;shiftdim(mean(fr(s2d31,wavesel,1:44),1),1)];
        s1_FR_6_t1=[s1_FR_6_t1;shiftdim(mean(fr(s1d61,wavesel,1:44),1),1)];
        s2_FR_6_t1=[s2_FR_6_t1;shiftdim(mean(fr(s2d61,wavesel,1:44),1),1)];
    end
    waveids=[waveids;waveidsess(wavesel)];

end


if opt.plot_1st_trial
    normfr=normalize([s1_FR_3,s2_FR_3,s1_FR_6,s2_FR_6,s1_FR_3_t1,s2_FR_3_t1,s1_FR_6_t1,s2_FR_6_t1].','range').';
    [s1_FR_3_t1,s2_FR_3_t1,s1_FR_6_t1,s2_FR_6_t1]=deal(normfr(:,(1:44)+44*4),normfr(:,(1:44)+44*5),normfr(:,(1:44)+44*6),normfr(:,(1:44)+44*7));
else
    normfr=normalize([s1_FR_3,s2_FR_3,s1_FR_6,s2_FR_6].','range').';
end
[s1_FR_3,s2_FR_3,s1_FR_6,s2_FR_6]=deal(normfr(:,1:44),normfr(:,(1:44)+44),normfr(:,(1:44)+88),normfr(:,(1:44)+132));

% save('FR_data.mat','s1_FR_3','s2_FR_3','s1_FR_6','s2_FR_6','s1_FR_3','s2_FR_3','s1_FR_6','s2_FR_6','s1_both',"s1_3only","s1_6only",'s2_both',"s2_3only","s2_6only")

% [coef,score,latent]=pca(normfr.');
% 
% [s1pc3,s2pc3,s1pc6,s2pc6]=deal(score(1:44,:),score((1:44)+44,:),score((1:44)+88,:),score((1:44)+132,:));

if opt.alt_proj
    %3s vs 6s axis
else
    cdMat=[s2_FR_3,s2_FR_6]-[s1_FR_3,s1_FR_6];
    cdDelay=mean(cdMat(:,[17:28,(17:40)+44]),2);
    cdDelay=cdDelay/norm(cdDelay);
    if opt.cd_pc_proj
        [qq,rr]=qr(cdDelay);
        cdfree=normfr.'*qq(:,2:end);
        [coef,score,latent]=pca(cdfree);
        [s1pc3,s2pc3,s1pc6,s2pc6]=deal(score(1:44,:),score((1:44)+44,:),score((1:44)+88,:),score((1:44)+132,:));
%     elseif opt.cd_tcom_proj
%         [s1pc3,s2pc3,s1pc6,s2pc6]=deal(s1_FR_3.'*comall,s2_FR_3.'*comall,s1_FR_6.'*comall,s2_FR_6.'*comall);
    end

%     [cdDelay3,cdDelay6,cdDelayB]=deal(cdDelay);
%     cdDelay3(~ismember(waveids,1:2))=0;
%     cdDelay6(~ismember(waveids,3:4))=0;
%     cdDelayB(~ismember(waveids,5:6))=0;
    
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

    
    if opt.cd_pc_proj
        [proj1X3,proj1Y3,proj1Z3]=deal((s1pc3(:,1)),(s1pc3(:,2)),(s1_FR_3.'*cdDelay));
        [proj2X3,proj2Y3,proj2Z3]=deal((s2pc3(:,1)),(s2pc3(:,2)),(s2_FR_3.'*cdDelay));
        [proj1X6,proj1Y6,proj1Z6]=deal((s1pc6(:,1)),(s2pc6(:,2)),(s1_FR_6.'*cdDelay));
        [proj2X6,proj2Y6,proj2Z6]=deal((s2pc6(:,1)),(s2pc6(:,2)),(s2_FR_6.'*cdDelay));
    elseif opt.cd_tcom_proj
        [proj1X3,proj1Y3,proj1Z3]=deal(s1_FR_3.'*comalls1./sum(s1_FR_3).',s1_FR_3.'*comalls2./sum(s1_FR_3).',(s1_FR_3.'*cdDelay));
        [proj2X3,proj2Y3,proj2Z3]=deal(s2_FR_3.'*comalls1./sum(s2_FR_3).',s2_FR_3.'*comalls2./sum(s2_FR_3).',(s2_FR_3.'*cdDelay));
        [proj1X6,proj1Y6,proj1Z6]=deal(s1_FR_6.'*comalls1./sum(s1_FR_6).',s1_FR_6.'*comalls2./sum(s1_FR_6).',(s1_FR_6.'*cdDelay));
        [proj2X6,proj2Y6,proj2Z6]=deal(s2_FR_6.'*comalls1./sum(s2_FR_6).',s2_FR_6.'*comalls2./sum(s2_FR_6).',(s2_FR_6.'*cdDelay));
    else
        [proj1X3,proj1Y3,proj1Z3]=deal(smooth(s1_FR_3.'*cdDelay3),smooth(s1_FR_3.'*cdDelay6),smooth(s1_FR_3.'*cdDelayB));
        [proj2X3,proj2Y3,proj2Z3]=deal(smooth(s2_FR_3.'*cdDelay3),smooth(s2_FR_3.'*cdDelay6),smooth(s2_FR_3.'*cdDelayB));
        [proj1X6,proj1Y6,proj1Z6]=deal(smooth(s1_FR_6.'*cdDelay3),smooth(s1_FR_6.'*cdDelay6),smooth(s1_FR_6.'*cdDelayB));
        [proj2X6,proj2Y6,proj2Z6]=deal(smooth(s2_FR_6.'*cdDelay3),smooth(s2_FR_6.'*cdDelay6),smooth(s2_FR_6.'*cdDelayB));
    end
    if opt.data_only
        data_str.proj1X3=proj1X3;
        data_str.proj1Y3=proj1Y3;
        data_str.proj1Z3=proj1Z3;
        data_str.proj1X6=proj1X6;
        data_str.proj1Y6=proj1Y6;
        data_str.proj1Z6=proj1Z6;
        data_str.proj2X3=proj2X3;
        data_str.proj2Y3=proj2Y3;
        data_str.proj2Z3=proj2Z3;
        data_str.proj2X6=proj2X6;
        data_str.proj2Y6=proj2Y6;
        data_str.proj2Z6=proj2Z6;
        return
    end

    if opt.plot_1st_trial
        [proj1E31,proj1M31,proj1L31]=deal(smooth(s1_FR_3_t1.'*cdDelay3),smooth(s1_FR_3_t1.'*cdDelay6),smooth(s1_FR_3_t1.'*cdDelayB));
        [proj2E31,proj2M31,proj2L31]=deal(smooth(s2_FR_3_t1.'*cdDelay3),smooth(s2_FR_3_t1.'*cdDelay6),smooth(s2_FR_3_t1.'*cdDelayB));
        [proj1E61,proj1M61,proj1L61]=deal(smooth(s1_FR_6_t1.'*cdDelay3),smooth(s1_FR_6_t1.'*cdDelay6),smooth(s1_FR_6_t1.'*cdDelayB));
        [proj2E61,proj2M61,proj2L61]=deal(smooth(s2_FR_6_t1.'*cdDelay3),smooth(s2_FR_6_t1.'*cdDelay6),smooth(s2_FR_6_t1.'*cdDelayB));
    end
    keyboard()
    if opt.gen_movie
        v=VideoWriter(sprintf('cd_pc_devp_movie_%s.mp4',opt.mem_type),'MPEG-4');
        open(v);
    end
    fh=figure('Color','w','Position',[100,100,900,600]);
    hold on;
    if strcmp(opt.mem_type,'3s')
        span=32;
    else
        span=44;
    end
       
    for tt=1:span
        cla
        switch opt.mem_type
            case '3s'
                view([5,-5])
            case '6s'
                view([6.5,10.5])
            case 'both'
                view([85,6.25])
            case 'mem'
                view([246,-5])
        end
        h13=plot3(proj1X3(1:tt),proj1Y3(1:tt),proj1Z3(1:tt),'-c','LineWidth',1.5);
        h23=plot3(proj2X3(1:tt),proj2Y3(1:tt),proj2Z3(1:tt),'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
        h16=plot3(proj1X6(1:tt),proj1Y6(1:tt),proj1Z6(1:tt),'-k','LineWidth',1.5);
        h26=plot3(proj2X6(1:tt),proj2Y6(1:tt),proj2Z6(1:tt),'-r','LineWidth',1.5);
        legend([h26,h23,h16,h13],{'S1-6s trial','S1-3s trial','S2-6s trial','S2-3s trial'},'FontSize',16,'Location','eastoutside','Orientation','vertical','AutoUpdate','off')

        th1=text(proj1X3(tt),proj1Y3(tt),proj1Z3(tt),'S2-3s','Color','c','FontSize',16);
        th2=text(proj2X3(tt),proj2Y3(tt),proj2Z3(tt),'S1-3s','Color',[0.9290 0.6940 0.1250],'FontSize',16);
        th3=text(proj1X6(tt),proj1Y6(tt),proj1Z6(tt),'S2-6s','Color','k','FontSize',16);
        th4=text(proj2X6(tt),proj2Y6(tt),proj2Z6(tt),'S1-6s','Color','r','FontSize',16);
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
        title(sprintf('%s selective, %s (%d msec)',opt.mem_type,tag,tt*250-4000),'FontSize',20)

        if opt.cd_pc_proj
            xlabel('PC1')
            ylabel('PC2')
            zlabel('CD')
%             xlim([-7,18])
%             ylim([-6,8])
%             zlim([-4,6])
        else
            xlabel('3s selective')
            ylabel('6s selective')
            zlabel('Both selective')
        end
        if opt.gen_movie
            for frame=1:15
                writeVideo(v,getframe(fh))
            end
        end
    end
    if opt.gen_movie
        close(v);
    end
    savefig(sprintf('cd_pc_traj_%s',opt.mem_type));

    % speed vs distance
    for tt=1:span
        dist3(tt)=norm(s1_FR_3(:,tt)-s2_FR_3(:,tt));
        dist6(tt)=norm(s1_FR_6(:,tt)-s2_FR_6(:,tt));
%         speed31(tt-1)=norm([proj1X3(tt),proj1Y3(tt),proj1Z3(tt)]-[proj1X3(tt-1),proj1Y3(tt-1),proj1Z3(tt-1)]);
% %         speed32
%         speed61(tt-1)=norm([proj1X6(tt),proj1Y6(tt),proj1Z6(tt)]-[proj1X6(tt-1),proj1Y6(tt-1),proj1Z6(tt-1)]);
%         speed62
    end

    %% dist vs speed
    if true
        fsh=figure('Color','w','Position',[32,32,800,400]);
        hold on
        d3h=plot(dist3,'-b','LineWidth',1);
        d6h=plot(dist6,'-r','LineWidth',1);
%         s3h=plot(speed31,'--b','LineWidth',1);
%         s6h=plot(speed61,'-b','LineWidth',1);
        set(gca,'XTick',3.5:12:span,'XTickLabel',-3:3:((span-20)/4))
        arrayfun(@(x) xline(x,':k'),[11.5,15.5,27.5,31.5,39.5,43.5])
        xlim([1,span]);
        xlabel('Time (s)')
        ylabel('S1-S2 euclidean distance (A.U.)')
        legend([d3h,d6h],{'S1/S2 distance, 3s Delay','S1/S2 distance, 6s Delay'},'Location','eastoutside')
        exportgraphics(fsh,sprintf('Dist_speed_%s.pdf',opt.mem_type),'ContentType','vector')
    end
    %%
    keyboard();

%     h13=plot3(proj1X3(1:tt),proj1Y3(1:tt),proj1Z3(1:tt),'-c','LineWidth',1.5);
%     h23=plot3(proj2X3(1:tt),proj2Y3(1:tt),proj2Z3(1:tt),'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
%     h16=plot3(proj1X6(1:tt),proj1Y6(1:tt),proj1Z6(1:tt),'-k','LineWidth',1.5);
%     h26=plot3(proj2X6(1:tt),proj2Y6(1:tt),proj2Z6(1:tt),'-r','LineWidth',1.5);

    delete(th1)
    delete(th2)
    delete(th3)
    delete(th4)

    th1=text(proj1X3(20),proj1Y3(20),proj1Z3(20),'S2-3s','Color','c','FontSize',16);
    th2=text(proj2X3(24),proj2Y3(24),proj2Z3(24),'S1-3s','Color',[0.9290 0.6940 0.1250],'FontSize',16);
    th3=text(proj1X6(28),proj1Y6(28),proj1Z6(28),'S2-6s','Color','k','FontSize',16);
    th4=text(proj2X6(32),proj2Y6(32),proj2Z6(32),'S1-6s','Color','r','FontSize',16);

    v=VideoWriter('cd_rotate_movie.mp4','MPEG-4');
    open(v);
    for dg=-5:10:355
        view([dg,10])
        for fm=1:10
            writeVideo(v,getframe(fh))
        end
    end
    close(v)

    if opt.plot_1st_trial
        h131=plot3(proj1E31(1:tt),proj1M31(1:tt),proj1L31(1:tt),'--c','LineWidth',1.5);
        h231=plot3(proj2E31(1:tt),proj2M31(1:tt),proj2L31(1:tt),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
        h161=plot3(proj1E61(1:tt),proj1M61(1:tt),proj1L61(1:tt),'--k','LineWidth',1.5);
        h261=plot3(proj2E61(1:tt),proj2M61(1:tt),proj2L61(1:tt),'--r','LineWidth',1.5);
    end
    v=VideoWriter('cd_PC_rotate_movie_blockN1.mp4','MPEG-4');
    open(v);
    for dg=120:10:480
        view([dg,30])
        for fm=1:10
            writeVideo(v,getframe(fh))
        end
    end
    close(v)

    distmat=[proj1X6-proj1X3,proj1Y6-proj1Y3,proj1Z6-proj1Z3];
    distmat1=[proj1E61-proj1E31,proj1M61-proj1M31,proj1L61-proj1L31];
    dist_3_6=arrayfun(@(x) norm(distmat(x,:)),1:44);
    dist_31_61=arrayfun(@(x) norm(distmat1(x,:)),1:44);

    distmat2=[proj2X6-proj2X3,proj2Y6-proj2Y3,proj2Z6-proj2Z3];
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


    distmat=[proj1X6-proj2X6,proj1Y6-proj2Y6,proj1Z6-proj2Z6];
    distmat1=[proj1E61-proj2E61,proj1M61-proj2M61,proj1L61-proj2L61];
    dist_3_6=arrayfun(@(x) norm(distmat(x,:)),1:44);
    dist_31_61=arrayfun(@(x) norm(distmat1(x,:)),1:44);

    distmat2=[proj2X6-proj2X3,proj2Y6-proj2Y3,proj2Z6-proj2Z3];
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
