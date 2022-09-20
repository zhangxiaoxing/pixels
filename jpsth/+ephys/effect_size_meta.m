%TODO brain region filter, olfaction filter.
function out=effect_size_meta()

%% gen data
curr_sess=-1;
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
homedir=ephys.util.getHomedir('type','raw');
out=struct();
out.cohen_d_olf=nan(numel(meta.allcid),3);
out.cohen_d_dur=nan(numel(meta.allcid),3);
for ii=1:numel(meta.allcid)
    if meta.sess(ii)~=curr_sess
        disp(meta.sess(ii))
        fpath=fullfile(homedir,ephys.sessid2path(meta.sess(ii)),"FR_All_1000.hdf5");
        fr=h5read(fpath,'/FR_All'); %Trial x SU x time-bin
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        fr_t_align=fr(:,:,5:7); % early/late delay
        curr_sess=meta.sess(ii);          
        csel=trials(:,9)~=0 & trials(:,10)~=0 & ismember(trials(:,8),[3 6]) & ismember(trials(:,5),[4 8]);
%         esel=trls(:,10)==0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        s1sel=csel & trials(:,5)==4;
        s2sel=csel & trials(:,5)==8;
        d3sel=csel & trials(:,8)==3;
        d6sel=csel & trials(:,8)==6;
    end
    suidx=(suid==meta.allcid(ii));
    % per-su effective size

    frs1=squeeze(fr_t_align(s1sel,suidx,:));
    frs2=squeeze(fr_t_align(s2sel,suidx,:));
    dmean=mean(frs1)-mean(frs2);
    poolstdfr=sqrt((var(frs1)+var(frs2))./2);
    out.cohen_d_olf(ii,:)=dmean./poolstdfr;

    frd3=squeeze(fr_t_align(d3sel,suidx,:));
    frd6=squeeze(fr_t_align(d6sel,suidx,:));
    dmean=mean(frd3)-mean(frd6);
    poolstdfr=sqrt((var(frd3)+var(frd6))./2);
    out.cohen_d_dur(ii,:)=dmean./poolstdfr;
end
end