function plot_std()
arguments

end
% persistent com_str onepath_ delay_ selidx_ decision_ rnd_half_ curve_
[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir();
[p_FR_3,np_FR_3,p_FR_6,np_FR_6]=deal([]);
figure('Color','w')
spid=1;
for mem_type=["Both","6s","3s"]
    for ii=reshape(cell2mat(sessmap.keys()),1,[])

        disp(ii)
        fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        wavesess=ones(size(suid))*double(ii);
        waveidsess=ephys.get_wave_id(wavesess,suid,'early',false,'ctx',true);

        s1d3=trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0;
        s2d3=trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0;
        s1d6=trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0;
        s2d6=trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0;

        switch lower(mem_type)
            case '3s'
                wavesel1=waveidsess==1;
                wavesel2=waveidsess==2;

            case '6s'
                wavesel1=waveidsess==3;
                wavesel2=waveidsess==4;

            case 'both'
                wavesel1=waveidsess==5;
                wavesel2=waveidsess==6;

        end

        if ~any(wavesel1|wavesel2)
            continue;
        end
%         fr=normalize(fr,2,'range');
            fr=fr+0.01;
          p_FR_3=[p_FR_3;shiftdim(std(fr(s1d3,wavesel1,1:44),1)./mean(fr(s1d3,wavesel1,1:44),1),1);shiftdim(std(fr(s2d3,wavesel2,1:44),1)./mean(fr(s1d3,wavesel2,1:44),1),1)];
        np_FR_3=[np_FR_3;shiftdim(std(fr(s2d3,wavesel1,1:44),1)./mean(fr(s2d3,wavesel1,1:44),1),1);shiftdim(std(fr(s1d3,wavesel2,1:44),1)./mean(fr(s2d3,wavesel2,1:44),1),1)];
          p_FR_6=[p_FR_6;shiftdim(std(fr(s1d6,wavesel1,1:44),1)./mean(fr(s1d6,wavesel1,1:44),1),1);shiftdim(std(fr(s2d6,wavesel2,1:44),1)./mean(fr(s1d6,wavesel2,1:44),1),1)];
        np_FR_6=[np_FR_6;shiftdim(std(fr(s2d6,wavesel1,1:44),1)./mean(fr(s2d6,wavesel1,1:44),1),1);shiftdim(std(fr(s1d6,wavesel2,1:44),1)./mean(fr(s2d6,wavesel2,1:44),1),1)];

%           p_FR_3=[p_FR_3;shiftdim(std(fr(s1d3,wavesel1,1:44),1),1);shiftdim(std(fr(s2d3,wavesel2,1:44),1),1)];
%         np_FR_3=[np_FR_3;shiftdim(std(fr(s2d3,wavesel1,1:44),1),1);shiftdim(std(fr(s1d3,wavesel2,1:44),1),1)];
%           p_FR_6=[p_FR_6;shiftdim(std(fr(s1d6,wavesel1,1:44),1),1);shiftdim(std(fr(s2d6,wavesel2,1:44),1),1)];
%         np_FR_6=[np_FR_6;shiftdim(std(fr(s2d6,wavesel1,1:44),1),1);shiftdim(std(fr(s1d6,wavesel2,1:44),1),1)];
        %         waveids=[waveids;waveidsess(wavesel1)];

    end
    subplot(3,2,spid)
    hold on
    psem=std(p_FR_3)./sqrt(size(p_FR_3,1));
    npsem=std(np_FR_3)./sqrt(size(np_FR_3,1));
    fill([1:44,44:-1:1],[mean(p_FR_3)+psem,fliplr(mean(p_FR_3)-psem)],[0.9290 0.6940 0.1250],'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none','FaceAlpha',0.15)
    fill([1:44,44:-1:1],[mean(np_FR_3)+npsem,fliplr(mean(np_FR_3)-npsem)],'c','FaceColor','c','EdgeColor','none','FaceAlpha',0.15)
    ph=plot(mean(p_FR_3),'-','Color',[0.9290 0.6940 0.1250]);
    nph=plot(mean(np_FR_3),'-c');
    set(gca,'XTick',4.5:12:32.5,'XTickLabel',-3:3:3)
    arrayfun(@(x) xline(x,':k'),[11.5,15.5,27.5,31.5]+1)
    xlim([0.5,32.5]);
    xlabel('Time (s)')
    ylabel('STD of Normalized FR')
    legend([ph,nph],{'Preferred trials','Nonpreferred trials'},'Location','northoutside','Orientation','horizontal')
    title(sprintf('%s selective, 3s WM delay',mem_type))
    spid=spid+1;
    subplot(3,2,spid)
    hold on
    psem=std(p_FR_6)./sqrt(size(p_FR_6,1));
    npsem=std(np_FR_6)./sqrt(size(np_FR_6,1));
    fill([1:44,44:-1:1],[mean(p_FR_6)+psem,fliplr(mean(p_FR_6)-psem)],'r','FaceColor','r','EdgeColor','none','FaceAlpha',0.15)
    fill([1:44,44:-1:1],[mean(np_FR_6)+npsem,fliplr(mean(np_FR_6)-npsem)],'b','FaceColor','b','EdgeColor','none','FaceAlpha',0.15)

    ph=plot(mean(p_FR_6),'-r');
    nph=plot(mean(np_FR_6),'-b');
    set(gca,'XTick',4.5:12:44.5,'XTickLabel',-3:3:6)
    arrayfun(@(x) xline(x,':k'),[11.5,15.5,39.5,43.5]+1)
    xlim([0.5,44.5]);
    xlabel('Time (s)')
    ylabel('STD of Normalized FR')
    title(sprintf('%s selective, 6s WM delay',mem_type))
    legend([ph,nph],{'Preferred trials','Nonpreferred trials'},'Location','northoutside','Orientation','horizontal')
    spid=spid+1;
end
sgtitle('Fano factor')
