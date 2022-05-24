function plot_cd_pref_np(opt)
arguments
    opt.plot_rpt (1,1) double {mustBePositive,mustBeInteger} =2
    opt.stats_rpt (1,1) double {mustBePositive,mustBeInteger} =5
    opt.plot_traj (1,1) logical = true
    opt.alt_proj (1,1) logical = false
    opt.cd_pc_proj (1,1) logical = true
    opt.plot_1st_trial (1,1) logical = false
    opt.gen_movie (1,1) logical = true
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'mem','3s','6s','both'})} = 'both'
end
% persistent com_str onepath_ delay_ selidx_ decision_ rnd_half_ curve_

[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir('type','raw');

[pref_FR_3,np_FR_3,pref_FR_6,np_FR_6,waveids]=deal([]);

for ii=reshape(cell2mat(sessmap.keys()),1,[])
    disp(ii)
    fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
    fr=h5read(fpath,'/FR_All');
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    dur_resp=behav.tag_block(trials,'wt',false);

    wavesess=ones(size(suid))*double(ii);
    waveidsess=ephys.get_wave_id(wavesess,suid,'early',false,'ctx',true);

    s1d3=trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0;
    s2d3=trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0;
    s1d6=trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0;
    s2d6=trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0;

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

    if ~any(waveselS1)
        continue;
    end

    pref_FR_3=[pref_FR_3;shiftdim(mean(fr(s1d3,waveselS1,1:44),1),1);shiftdim(mean(fr(s2d3,waveselS2,1:44),1),1)];
    np_FR_3=[np_FR_3;shiftdim(mean(fr(s2d3,waveselS1,1:44),1),1);shiftdim(mean(fr(s1d3,waveselS2,1:44),1),1)];
    pref_FR_6=[pref_FR_6;shiftdim(mean(fr(s1d6,waveselS1,1:44),1),1);shiftdim(mean(fr(s2d6,waveselS2,1:44),1),1)];
    np_FR_6=[np_FR_6;shiftdim(mean(fr(s2d6,waveselS1,1:44),1),1);shiftdim(mean(fr(s1d6,waveselS2,1:44),1),1)];

    waveids=[waveids;waveidsess(waveselS1)];
end

normfr=normalize([pref_FR_3,np_FR_3,pref_FR_6,np_FR_6].','range').';
[pref_FR_3,np_FR_3,pref_FR_6,np_FR_6]=deal(normfr(:,1:44),normfr(:,(1:44)+44),normfr(:,(1:44)+88),normfr(:,(1:44)+132));

vecdiff3=pref_FR_3-np_FR_3;
vecdiff6=pref_FR_6-np_FR_6;

fh=figure('Color','w','Position',[32,32,250,250]);
hold on
% norm
% h3=plot(arrayfun(@(x) norm(vecdiff3(:,x)-vecdiff3(:,x-1)),2:size(vecdiff3,2)),'-b');
% h6=plot(arrayfun(@(x) norm(vecdiff6(:,x)-vecdiff6(:,x-1)),2:size(vecdiff6,2)),'-r');

% raw sum
% h3=plot(arrayfun(@(x) sum(vecdiff3(:,x)-vecdiff3(:,x-1)),2:size(vecdiff3,2)),'-b');
% h6=plot(arrayfun(@(x) sum(vecdiff6(:,x)-vecdiff6(:,x-1)),2:size(vecdiff6,2)),'-r');

% vector angle
h3=plot(arrayfun(@(x) ND_angle(vecdiff3(:,x),vecdiff3(:,x-1)),2:size(vecdiff3,2)),'-b');
h6=plot(arrayfun(@(x) ND_angle(vecdiff6(:,x),vecdiff6(:,x-1)),2:size(vecdiff6,2)),'-r');

set(gca(),'XTick',[4,16,28,40]-0.5,'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,':k'),[12,16,28,32,40,44]-0.5)
xlabel('Time (s)')
ylabel('Vector angular speed (degree)')
legend([h3,h6],{'In 3s delay','In 6s delay'},'Location','northoutside','Orientation','horizontal')
sgtitle(sprintf('%s selective',opt.mem_type))
exportgraphics(fh,sprintf('vector_angular_speed_%s.pdf',opt.mem_type),'ContentType','vector')

end

function deg=ND_angle(u,v)
rad = acos((u'*v)/(norm(u)*norm(v)));
deg = rad/pi*180;
end