% function wf_stats=waveform_dist()
% Extracellular Spike Waveform Dissociates Four Functionally Distinct Cell Classes in Primate Cortex
% Current Biology 2019 Earl K. Miller,Markus Siegel
% excluded waveforms that satisfied any of three criteria for atypical shape:
% (1) the amplitude of the main trough was smaller than the subsequent positive
% peak (n = 41), (2) the trace was noisy, defined as > = 6 local maxima of magnitude > = 0.01 (n = 38), (3) there was one or more local
% maxima in the period between the main trough and the subsequent peak (n = 35).

% aligned to baseline during extraction

[~,~,sessmap]=ephys.sessid2path(0);
homedir=replace(ephys.util.getHomedir(),'SPKINFO','WF');

wf_stats=[];
for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
    folder=sessmap(sesskey);
    seppath=regexp(fullfile(homedir,folder),'(^.*)(imec[01])(.*$)','tokens');
    folder0=[seppath{1}{1},'imec0',seppath{1}{3}];
    folder1=[seppath{1}{1},'imec1',seppath{1}{3}];
    [cid0,cid1,judge0,judge1]=deal([]);
    if isfolder(folder0),stats0=oneFolder(folder0,sesskey);end
    if isfolder(folder1),stats1=oneFolder(folder1,sesskey);end
    wf_stats=[wf_stats;stats0;stats1];
end
save('wf_stats.mat','wf_stats');
% end
for ii=reshape(find(wf_stats(:,7)<50 & wf_stats(:,3)>0),1,[])
    fh=figure();plot(wf_stats(ii,11:101));waitfor(fh);
end


function stats=oneFolder(folder,sess)

if ~isfile(fullfile(folder,'waveform.mat'))
    disp('Missing waveform file');
    disp(folder);
    stats=[];
    keyboard();
    return
end
fstr=load(fullfile(folder,'waveform.mat'));
wf_all=fstr.waveform;
stats=nan(size(wf_all,1),101);
if contains(folder,'imec1') && max([wf_all{:,2}])<10000
    warning('Missing probe# tag');
    wf_all(:,2)=num2cell([wf_all{:,2}].'+10000);
end

for ii=1:size(wf_all,1)
    stats(ii,1)=sess;
    stats(ii,2)=wf_all{ii,2};%cid
    wf=wf_all{ii,4};
    stats(ii,11:101)=wf;
    %criteria 1
    [mmin,mi]=min(wf);
    [mmax,xi]=max(wf((mi+1):end));
    if mmax>-mmin
        stats(ii,3)=-1;% inverted shape
        wf=-wf;
    else
        stats(ii,3)=1;
    end
    %criteria 2
    [lc_pk,lc_pk_ts]=findpeaks(wf,'MinPeakProminence',-0.25*mmin);
    stats(ii,4)=numel(lc_pk); % noise,

    %criteria 3
    if ~isempty(xi) && xi>3
        [lc_pk,lc_pk_ts]=findpeaks(wf((mi+1):(mi+xi-1)),'MinPeakProminence',-0.05*mmin);
        stats(ii,5)=numel(lc_pk);%noise after anti peak
    elseif isempty(xi) || (xi<=3) || (~isempty(lc_pk))
        stats(ii,5)=-1; % missing character current peak, should reject.
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

    if numel(lcross)==1 && numel(rcross)==1
        fwhm=rcross-lcross;
        stats(ii,6:7)=[deltaTS,fwhm];
    else
        stats(ii,6:7)=[-1,-1];
    end
end
% keyboard()
% out=cell2mat(wf_all(stats(:,1)>0,2));
% judge=stats(:,1);
end

