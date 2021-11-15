function plot_cd(opt)
arguments
    opt.keep_sust (1,1) logical = false % use sustained coding neuron
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6 % DPA delay duration
    opt.plot_rpt (1,1) double {mustBePositive,mustBeInteger} =10
    opt.stats_rpt (1,1) double {mustBePositive,mustBeInteger} =250
    opt.to_plot (1,1) logical = false
end
% persistent com_str onepath_ delay_ selidx_ decision_ rnd_half_ curve_

if opt.delay==6
    warning('Delay set to default 6')
end

meta_str=ephys.util.load_meta('type','neupix','delay',opt.delay);
homedir=ephys.util.getHomedir('type','raw');
fl=dir(fullfile(homedir,'**','FR_All_ 250.hdf5'));
[s1mat,s2mat]=deal([]);
[s1shufmat,s2shufmat]=deal(cell(1,2*opt.stats_rpt));

for ii=1:size(fl,1)
    dpath=regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once');
    fpath=fullfile(fl(ii).folder,fl(ii).name);
    fpath=replace(fpath,'\',filesep());
    pc_stem=replace(dpath,'/','\');
    sesssel=startsWith(meta_str.allpath,pc_stem);
    if ~any(sesssel), continue;end
    fr=h5read(fpath,'/FR_All');
    trial=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    %TODO nonmem,incongruent
    if opt.keep_sust
        mcid=meta_str.allcid(meta_str.mem_type>0 & sesssel.' & strcmp(meta_str.reg_tree(2,:),'CTX'));
    else
        mcid=meta_str.allcid((meta_str.mem_type==2 | meta_str.mem_type==4) & sesssel.' & strcmp(meta_str.reg_tree(2,:),'CTX'));
    end
    if isempty(mcid), continue; end
    msel=ismember(suid,mcid);
    s1sel=trial(:,5)==4 & trial(:,8)==opt.delay & trial(:,9)>0 & trial(:,10)>0;
    s2sel=trial(:,5)==8 & trial(:,8)==opt.delay & trial(:,9)>0 & trial(:,10)>0;

    for rpt=1:opt.stats_rpt
        s1pool=find(s1sel);
        s2pool=find(s2sel);
        s1a=randsample(s1pool,floor(numel(s1pool)./2));
        s1b=s1pool(~ismember(s1pool,s1a));
        s2a=randsample(s2pool,floor(numel(s2pool)./2));
        s2b=s2pool(~ismember(s2pool,s2a));
        s1shufmat{rpt*2-1}=[s1shufmat{rpt*2-1};shiftdim(mean(fr(s1a,msel,9:(opt.delay*4+16))),1)];
        s1shufmat{rpt*2}=[s1shufmat{rpt*2};shiftdim(mean(fr(s1b,msel,9:(opt.delay*4+16))),1)];
        s2shufmat{rpt*2-1}=[s2shufmat{rpt*2-1};shiftdim(mean(fr(s2a,msel,9:(opt.delay*4+16))),1)];
        s2shufmat{rpt*2}=[s2shufmat{rpt*2};shiftdim(mean(fr(s2b,msel,9:(opt.delay*4+16))),1)];
    end
    s1mat=[s1mat;shiftdim(mean(fr(s1sel,msel,9:(opt.delay*4+16))),1)];
    s2mat=[s2mat;shiftdim(mean(fr(s2sel,msel,9:(opt.delay*4+16))),1)];
end

cdMat=s1mat-s2mat;
cdDelayE=mean(cdMat(:,9:12),2);
cdDelayE=cdDelayE/norm(cdDelayE);

cdDelayM=mean(cdMat(:,19:22),2);
cdDelayM=cdDelayM/norm(cdDelayM);

cdDelayL=mean(cdMat(:,29:32),2);
cdDelayL=cdDelayL/norm(cdDelayL);
if opt.to_plot
    fh=figure('Color','w','Position',[100,100,190,190]);
    hold on;
    for rpt=1:(opt.plot_rpt*2)
        [proj1E,proj1M,proj1L]=deal(s1shufmat{rpt}.'*cdDelayE,s1shufmat{rpt}.'*cdDelayM,s1shufmat{rpt}.'*cdDelayL);
        [proj2E,proj2M,proj2L]=deal(s2shufmat{rpt}.'*cdDelayE,s2shufmat{rpt}.'*cdDelayM,s2shufmat{rpt}.'*cdDelayL);
        %     plot3(proj1E(1:5),proj1M(1:5),proj1L(1:5),':r.','LineWidth',0.5)
        %     plot3(proj2E(1:5),proj2M(1:5),proj2L(1:5),':b.','LineWidth',0.5)
        %     plot3(proj1E(5:9),proj1M(5:9),proj1L(5:9),'--r.','LineWidth',0.5)
        %     plot3(proj2E(5:9),proj2M(5:9),proj2L(5:9),'--b.','LineWidth',0.5)
        plot3(smooth(proj1E,3),smooth(proj1M,3),smooth(proj1L,3),'-m','LineWidth',0.5)
        plot3(smooth(proj2E,3),smooth(proj2M,3),smooth(proj2L,3),'-c','LineWidth',0.5)
        %     plot3(proj1E(9),proj1M(9),proj1L(9),'ko','MarkerFaceColor','k')
        %     plot3(proj2E(9),proj2M(9),proj2L(9),'ko','MarkerFaceColor','k')
    end
    % plot3(proj1E(1:5),proj1M(1:5),proj1L(1:5),':r.','LineWidth',1)
    % plot3(proj2E(1:5),proj2M(1:5),proj2L(1:5),':b.','LineWidth',1)

    [proj1E,proj1M,proj1L]=deal(smooth(s1mat.'*cdDelayE),smooth(s1mat.'*cdDelayM),smooth(s1mat.'*cdDelayL));
    [proj2E,proj2M,proj2L]=deal(smooth(s2mat.'*cdDelayE),smooth(s2mat.'*cdDelayM),smooth(s2mat.'*cdDelayL));


    plot3(proj1E(1:9),proj1M(1:9),proj1L(1:9),':r.','LineWidth',1)
    plot3(proj2E(1:9),proj2M(1:9),proj2L(1:9),':b.','LineWidth',1)
    plot3(proj1E(9:end),proj1M(9:end),proj1L(9:end),'-r.','LineWidth',1)
    plot3(proj2E(9:end),proj2M(9:end),proj2L(9:end),'-b.','LineWidth',1)

    plot3(proj1E(9),proj1M(9),proj1L(9),'ko','MarkerFaceColor','k')
    plot3(proj2E(9),proj2M(9),proj2L(9),'ko','MarkerFaceColor','k')

    xlabel('Early CD (A.U.)')
    ylabel('Mid CD (A.U.)')
    zlabel('Late CD (A.U.)')
    grid on
    keyboard()
    exportgraphics(fh,'CD_project_traj_3d.pdf','ContentType','vector');
end
estat=[];
mstat=[];
lstat=[];

for rpt=1:opt.stats_rpt
    [proj1E,proj1M,proj1L]=deal(s1shufmat{rpt}.'*cdDelayE,s1shufmat{rpt}.'*cdDelayM,s1shufmat{rpt}.'*cdDelayL);
    [proj2E,proj2M,proj2L]=deal(s2shufmat{rpt}.'*cdDelayE,s2shufmat{rpt}.'*cdDelayM,s2shufmat{rpt}.'*cdDelayL);
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
