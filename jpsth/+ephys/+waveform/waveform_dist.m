close all
figure('Color','w');
stats=[];% type, trough-peak, fwhm
showcaseCount=0;
fl=dir('K:\neupix\WF\neuropixel\**\waveform.mat');
for f=fl'
    fstr=load(fullfile(f.folder,f.name));
    wf_all=fstr.waveform;
    for one_wf=wf_all(:,4)'
        wf=one_wf{1};
        %criteria 1
        if max(wf)>-min(wf)
            disp('type 1 bad wf')
            stats(end+1,:)=[-1,0,0];
            continue
        end
        %criteria 2
        [lc_pk,lc_pk_ts]=findpeaks(wf,'MinPeakProminence',-0.05*min(wf));
        if numel(lc_pk)>=6
            disp('type 2 bad wf')
%             findpeaks(wf,'MinPeakProminence',-0.05*min(wf));
            stats(end+1,:)=[-2,0,0];
            continue
        end
        %criteria 3
        [~,t_ts]=min(wf);
        [~,p_ts]=max(wf(t_ts+1:end)); 
        if p_ts>3
            [lc_pk,lc_pk_ts]=findpeaks(wf(t_ts:(t_ts+p_ts-1)),'MinPeakProminence',-0.05*min(wf)); 
        end
        if isempty(p_ts) || (p_ts<=3) || (~isempty(lc_pk))
            disp('type 3 bad wf')
            if isempty(p_ts) || (p_ts<=3)
%                 plot(wf)
            else
%                 findpeaks(wf(t_ts:(t_ts+p_ts-1)),'MinPeakProminence',-0.05*min(wf))
            end
            stats(end+1,:)=[-3,0,0];
            pause(0.2)
            continue
        end
        %trough_peak dist
        wf=spline(1:91,wf,1:0.03:91);
        scale=max(abs(wf));
        wf=wf./scale;
        [trough,troughTS]=min(wf);
        [lp,deltaTS]=max(wf((troughTS+1):end));%to late peak
        
        %fwhm
        lcross=find(wf<-0.5,1);
        rcross=find(wf(troughTS:end)>-0.5,1)+troughTS;
%         figure()
        if numel(lcross)==1 && numel(rcross)==1
            fwhm=rcross-lcross;
            stats(end+1,:)=[0,deltaTS,fwhm];
            showcaseCount=showcaseCount+1;

            if showcaseCount>200 && showcaseCount<=250
                hold on
                plot(wf)
                plot(lcross,-0.5,'ro')
                plot(rcross,-0.5,'ro')
                plot(troughTS,trough,'bo')
                plot(troughTS+deltaTS,lp,'bo')
                hold off
            end
        else
            disp('type 4 bad wf')
%             plot(wf)
            stats(end+1,:)=[-4,0,0];
            pause(0.2)
            continue
        end
    end
end

% scatter(stats(sel,2),stats(sel,3),5,'k','o','MarkerFaceColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none')
% xlim([0,1250])
% ylim([0,600])
histogram(stats(sel,2),0:20:1100)
% histogram(stats(sel,3),50:10:400)
fh=gcf();
fh.Color='w';
set(gca(),'YTick',0:500:2000,'XTick',0:500:1000)
xlabel('trough to peak (us)')
ylabel('number of S.U.')
print('waveform_hist.eps','-depsc','-painters')
