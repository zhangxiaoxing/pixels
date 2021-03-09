function cid=goodCid(folder,opt)
arguments
    folder (1,:) char
    opt.ntrack (1,1) double {mustBeMember(opt.ntrack,[1,2])} = 2
    opt.waveform (1,1) logical = true
end
cid=oneFolder(folder);
if opt.ntrack>1
    cid1=oneFolder(replace(folder,'imec0','imec1'));
    cid=[cid;cid1+10000];
end
end


function cluster_ids=oneFolder(folder)
s1s=30000;
FR_Th=1.0;
metaf=dir(fullfile(folder,'*.meta'));
fh=fopen(fullfile(metaf.folder,metaf.name));
ts=textscan(fh,'%s','Delimiter',{'\n'});
nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
spkNThresh=nSample/385/s1s/2*FR_Th;
clusterInfo = readtable(fullfile(folder,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
contamGood=strcmp(clusterInfo{:,4},'good');
freqGood=clusterInfo{:,10}>spkNThresh;
cluster_ids = table2array(clusterInfo(contamGood & freqGood,1));

end