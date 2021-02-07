function [path0,path1,error_list]=jointFolder(folder,error_list,homedir)
    sessFolder=replace(fullfile(fullfile(homedir,'DataSum'),folder),'\',filesep);
    spkpath0=dir(fullfile(sessFolder,'*imec0*','spike_info.mat'));
    if isempty(spkpath0)
        sessFolder=replace(sessFolder,'DataSum','DataSum/singleProbe');
        spkpath0=dir(fullfile(sessFolder,'*imec0*','spike_info.mat'));
    end
    if isempty(spkpath0)
        path0=[];
    else
        path0=spkpath0.folder;
    end
    spkpath1=dir(fullfile(sessFolder,'*imec1*','spike_info.mat'));
    if isempty(spkpath1)
        path1=[];
    else
        path1=spkpath1.folder;
    end
end



