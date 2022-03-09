fl=dir('K:\neupix\META\**\cluster_info.tsv');
amps=[];
for ii=1:numel(fl)
    disp(ii)
    tbl=readtable(fullfile(fl(ii).folder,fl(ii).name),'FileType','text');
    if any(strcmp(tbl.Properties.VariableNames,'fr'))
        fr=tbl.fr>1.0;
    else
        fr=tbl.firing_rate>1.0;
    end
    contam=strcmp(tbl.KSLabel,'good');
    amp=tbl.Amplitude(fr&contam);
    amps=[amps;amp];
end