function dist_to_attactor(opt)
arguments
    opt.plot_rpt (1,1) double {mustBePositive,mustBeInteger} =2
    opt.stats_rpt (1,1) double {mustBePositive,mustBeInteger} =5
    opt.plot_traj (1,1) logical = true
    opt.alt_proj (1,1) logical = false
    opt.cd_pc_proj (1,1) logical = true
    opt.plot_1st_trial (1,1) logical = false
    opt.gen_movie (1,1) logical = true
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'sust','3s','6s','both'})} = 'sust'
end
% persistent com_str onepath_ delay_ selidx_ decision_ rnd_half_ curve_
meta6=ephys.util.load_meta('delay',6);
[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir();

[pref_FR_3,np_FR_3,pref_FR_6,np_FR_6]=deal([]);

for ii=reshape(cell2mat(sessmap.keys()),1,[])
    disp(ii)
    fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
    fr=h5read(fpath,'/FR_All');
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    s1d3=trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0;
    s2d3=trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0;
    s1d6=trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0;
    s2d6=trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0;

    if ~strcmp(opt.mem_type,'sust')
        wavesess=ones(size(suid))*double(ii);
        waveidsess=ephys.get_wave_id(wavesess,suid,'early',false,'ctx',true);

        switch opt.mem_type
            case 'mem'
                waveselS1=ismember(waveidsess,[1 3 5]);
                waveselS2=ismember(waveidsess,[2 4 6]);
            case '3s'
                waveselS1=waveidsess==1;
                waveselS2=waveidsess==2;
            case '6s'
                waveselS1=waveidsess==3;
                waveselS2=waveidsess==4;
            case 'both'
                waveselS1=waveidsess==5;
                waveselS2=waveidsess==6;
        end
    else
        waveselS1=meta6.mem_type(meta6.sess==ii)==1;
        waveselS2=meta6.mem_type(meta6.sess==ii)==3;
    end

    if ~any(waveselS1 | waveselS2)
        continue;
    end

    pref_FR_3=[pref_FR_3;shiftdim(mean(fr(s1d3,waveselS1,1:44),1),1);shiftdim(mean(fr(s2d3,waveselS2,1:44),1),1)];
    np_FR_3=[np_FR_3;shiftdim(mean(fr(s2d3,waveselS1,1:44),1),1);shiftdim(mean(fr(s1d3,waveselS2,1:44),1),1)];
    pref_FR_6=[pref_FR_6;shiftdim(mean(fr(s1d6,waveselS1,1:44),1),1);shiftdim(mean(fr(s2d6,waveselS2,1:44),1),1)];
    np_FR_6=[np_FR_6;shiftdim(mean(fr(s2d6,waveselS1,1:44),1),1);shiftdim(mean(fr(s1d6,waveselS2,1:44),1),1)];

%     waveids=[waveids;waveidsess(waveselS1)];
end

normfr=normalize([pref_FR_3,np_FR_3,pref_FR_6,np_FR_6].','range').';
[pref_FR_3,np_FR_3,pref_FR_6,np_FR_6]=deal(normfr(:,1:44),normfr(:,(1:44)+44),normfr(:,(1:44)+88),normfr(:,(1:44)+132));


fh=figure('Color','w','Position',[1,31,1440,793]);
if ~strcmp(opt.mem_type,'sust')
subplot(2,2,1)
hold on
for cdBin=1:3
    coding_attractor=mean(pref_FR_3(:,(14:15)+cdBin*4)-np_FR_3(:,(14:15)+cdBin*4),2);
    cd_dist=arrayfun(@(x) norm(pref_FR_3(:,x)-np_FR_3(:,x)-coding_attractor),1:size(np_FR_3,2));
    plot(cd_dist)
end
set(gca(),'XTick',[4,16,28,40]+0.5,'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,':k'),[12,16,28,32]+0.5)
xlabel('Time (s)')
ylabel('Distance to hypothetical attractor')
title('3s trial, distance')
subplot(2,2,2)
hold on
for cdBin=1:3
    coding_attractor=mean(pref_FR_3(:,(14:15)+cdBin*4)-np_FR_3(:,(14:15)+cdBin*4),2);
    cd_dist=arrayfun(@(x) norm(pref_FR_3(:,x)-np_FR_3(:,x)-coding_attractor),1:size(np_FR_3,2));
    plot(-diff(cd_dist))
end
set(gca(),'XTick',[4,16,28,40],'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,':k'),[12,16,28,32])
xlabel('Time (s)')
ylabel('Speed toward hypothetical attractor')
title('3s trial, speed')
end
subplot(2,2,3)
hold on
for cdBin=1:6
    coding_attractor=mean(pref_FR_6(:,(14:15)+cdBin*4)-np_FR_6(:,(14:15)+cdBin*4),2);
    cd_dist=arrayfun(@(x) norm(pref_FR_6(:,x)-np_FR_6(:,x)-coding_attractor),1:size(np_FR_6,2));
    plot(cd_dist)
end

set(gca(),'XTick',[4,16,28,40]+0.5,'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,':k'),[12,16,40,44]+0.5)
xlabel('Time (s)')
ylabel('Distance to hypothetical attractor')
title('6s trial, distance')

subplot(2,2,4)
hold on
for cdBin=1:6
    coding_attractor=mean(pref_FR_6(:,(14:15)+cdBin*4)-np_FR_6(:,(14:15)+cdBin*4),2);
    cd_dist=arrayfun(@(x) norm(pref_FR_6(:,x)-np_FR_6(:,x)-coding_attractor),1:size(np_FR_6,2));
    plot(-diff(cd_dist))
end
set(gca(),'XTick',[4,16,28,40],'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,':k'),[12,16,40,44])
xlabel('Time (s)')
ylabel('Speed toward hypothetical attractor')
title('6s trial, speed')
sgtitle(sprintf('Distance and speed to hypothetical attactors, %s selective',opt.mem_type))
exportgraphics(fh,sprintf('dist_speed_attractor_%s.pdf',opt.mem_type),'ContentType','vector')

end
