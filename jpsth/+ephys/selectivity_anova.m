meta=ephys.util.load_meta();
[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir('type','raw');
anovameta=struct();
[anovameta.sess,anovameta.allcid,anovameta.anovap]=deal([]);
for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
    disp(sesskey)
    fpath=fullfile(homedir,sessmap(sesskey),"FR_All_1000.hdf5");
    fr=h5read(fpath,'/FR_All');
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');
    
    sensible_sel=ismember(trials(:,8),[3,6]) & all(trials(:,9:10)~=0,2) & ismember(trials(:,5),[4,8]);
    
    trials=trials(sensible_sel,:);
    fr=fr(sensible_sel,:,:);
    
    anovanps=nan(size(fr,2),7); %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin

    for su=1:size(fr,2)
        frmat=squeeze(fr(:,su,4:7));
        sufrvec=frmat(:);
        sampvec=repmat(trials(:,5),4,1);
        durvec=repmat(trials(:,8),4,1);
        binvec=reshape(repmat((1:4),size(fr,1),1),[],1);
        anovanps(su,:)=anovan(sufrvec,{sampvec,durvec,binvec},'model','full','display','off');
    end

    anovameta.sess=[anovameta.sess;repmat(sesskey,size(fr,2),1)];
    anovameta.allcid=[anovameta.allcid;suid];
    anovameta.anovap=[anovameta.anovap;anovanps];
end
if false
    save('anovameta.mat','anovameta')
end
keyboard()
