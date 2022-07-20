function prefwave=pref_wave()
prefwave=struct();

su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
prefwave.pref_id=zeros(size(su_meta.allcid));
prefwave.fr_13_23_16_26=zeros(numel(su_meta.allcid),4);
currstem=[];
homedir=ephys.util.getHomedir('type','raw');
for ii=1:numel(su_meta.allcid)
    if ~strcmp(currstem,su_meta.allpath{ii})
        currstem=su_meta.allpath{ii};
        fpath=replace(fullfile(homedir,su_meta.allpath{ii},'FR_All_1000.hdf5'),'\',filesep());
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        fr=h5read(fpath,'/FR_All');
%         disp(num2str(ii)+" total trials "+num2str(nnz(trials(:,9) & trials(:,10))))
        s1d3t=(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==3);
        s2d3t=(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==3);
        s1d6t=(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6);
        s2d6t=(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6);
    end
    
    
    mm=[mean(fr(s1d3t,suid==su_meta.allcid(ii),5:7),'all'),...
        mean(fr(s2d3t,suid==su_meta.allcid(ii),5:7),'all'),...
        mean(fr(s1d6t,suid==su_meta.allcid(ii),5:7),'all'),...
        mean(fr(s2d6t,suid==su_meta.allcid(ii),5:7),'all')];
    [~,mmidx]=max(mm);
    prefwave.pref_id(ii)=mmidx;
    prefwave.fr_13_23_16_26(ii,:)=mm;
end
end