function cid=goodWaveform(folder,opt)

% Extracellular Spike Waveform Dissociates Four Functionally Distinct Cell Classes in Primate Cortex
% Current Biology 2019 Earl K. Miller,Markus Siegel
% excluded waveforms that satisfied any of three criteria for atypical shape:
% (1) the amplitude of the main trough was smaller than the subsequent positive
% peak (n = 41), (2) the trace was noisy, defined as > = 6 local maxima of magnitude > = 0.01 (n = 38), (3) there was one or more local
% maxima in the period between the main trough and the subsequent peak (n = 35).

% aligned to baseline during extraction


arguments
    folder (1,:) char
    opt.presel (1,:) double = []
    opt.stats (1,1) logical = false
end
seppath=regexp(folder,'(^.*)(imec[01])(.*$)','tokens');
folder0=[seppath{1}{1},'imec0',seppath{1}{3}];
folder1=[seppath{1}{1},'imec1',seppath{1}{3}];
cid0=[];cid1=[];
if isfolder(folder0),cid0=oneFolder(folder0,opt);end
if isfolder(folder1),cid1=oneFolder(folder1,opt);end
cid=[cid0;cid1+10000];
end


function out=oneFolder(folder,opt)
stats=[];% type, trough-peak, fwhm
showcaseCount=0;
if ~isfile(fullfile(folder,'waveform.mat'))
    disp('Missing waveform file');
    disp(folder);
    out=[];
    keyboard();
    return
end
fstr=load(fullfile(folder,'waveform.mat'));
wf_all=fstr.waveform;
% for one_wf=wf_all(:,4)'
for i=1:size(wf_all,1)
    cid=wf_all{i,2};
    if ~isempty(opt.presel) && ~ismember(cid,opt.presel)
        stats(end+1,:)=[-5,0,0];
        continue
    end
    wf=wf_all{i,4};
    %criteria 1
    [mmin,mi]=min(wf);
    [mmax,xi]=max(wf((mi+1):end));
    if mmax>-mmin
%         disp('type 1 bad wf')
        stats(end+1,:)=[-1,0,0];
        continue
    end
    %criteria 2
    [lc_pk,lc_pk_ts]=findpeaks(wf,'MinPeakProminence',-0.05*mmin);
    if numel(lc_pk)>=6
%         disp('type 2 bad wf')
        stats(end+1,:)=[-2,0,0];
        continue
    end
    %criteria 3
    if ~isempty(xi) && xi>3
        [lc_pk,lc_pk_ts]=findpeaks(wf((mi+1):(mi+xi-1)),'MinPeakProminence',-0.05*mmin);
    end
    if isempty(xi) || (xi<=3) || (~isempty(lc_pk))
%         disp('type 3 bad wf')
        stats(end+1,:)=[-3,0,0];
        continue
    end
    stats(end+1,:)=[1,0,0];
    if opt.stats
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
out=cell2mat(wf_all(stats(:,1)>0,2));

end
