function cid=goodCid(folder,opt)
arguments
    folder (1,:) char
    opt.FR_th (1,1) double = 1.0
end
seppath=regexp(folder,'(^.*)(imec[01])(.*$)','tokens');
folder0=[seppath{1}{1},'imec0',seppath{1}{3}];
folder1=[seppath{1}{1},'imec1',seppath{1}{3}];
cid0=[];cid1=[];
if isfolder(folder0),cid0=oneFolder(folder0,opt);end
if isfolder(folder1),cid1=oneFolder(folder1,opt);end
cid=[cid0;cid1+10000];
end


function cluster_ids=oneFolder(folder,opt)
s1s=30000;
metaf=dir(fullfile(folder,'*.meta'));
fh=fopen(fullfile(metaf.folder,metaf.name));
ts=textscan(fh,'%s','Delimiter',{'\n'});
nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
spkNThresh=nSample/385/s1s/2*opt.FR_th;
clusterInfo = readtable(fullfile(folder,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
contamGood=strcmp(clusterInfo{:,4},'good');
freqGood=clusterInfo{:,10}>spkNThresh;
cluster_ids = table2array(clusterInfo(contamGood & freqGood,1));

end