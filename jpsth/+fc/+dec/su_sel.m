function out=su_sel(sums,folder,delay,opt)
arguments
    sums (:,3) cell
    folder (1,:) char
    delay (1,1) double
    opt.type (1,:) char {mustBeMember(opt.type,{'all','pre','post','none','both'})} = 'all'
end

if strcmp(opt.type,'all')
    out=sums;
    return
else
    sus_trans=h5read(sprintf('../transient_%d.hdf5',delay),'/sus_trans');%4+2+7
    allpath=h5read(sprintf('../transient_%d.hdf5',delay),'/path');
    allcid=h5read(sprintf('../transient_%d.hdf5',delay),'/cluster_id');
end
% keyboard()

fstub=regexp(folder,'(?<=.*DataSum\\)[^\\]*','match','once');
sel0=startsWith(allpath,fstub) & contains(allpath,'imec0');
sel1=startsWith(allpath,fstub) & contains(allpath,'imec1');
pairs=cell2mat(sums(:,2)');

nonsel=[allcid(~any(sus_trans,2) & sel0);allcid(~any(sus_trans,2) & sel1)+10000];
sel=[allcid(any(sus_trans(:,[1,2,4]),2) & sel0);allcid(any(sus_trans(:,[1,2,4]),2) & sel1)+10000];    

if strcmp(opt.type,'none')
    susel=(ismember(pairs(1,:),nonsel) & ismember(pairs(2,:),nonsel))';
elseif strcmp(opt.type,'both')
    susel=(ismember(pairs(1,:),sel) & ismember(pairs(2,:),sel))';
elseif strcmp(opt.type,'pre')
    susel=(ismember(pairs(1,:),sel) & ismember(pairs(2,:),nonsel))';
elseif strcmp(opt.type,'post')
    susel=(ismember(pairs(1,:),nonsel) & ismember(pairs(2,:),sel))';
end

out=sums(susel,:);

end