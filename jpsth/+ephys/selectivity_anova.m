function anovameta_=selectivity_anova(opt)
arguments
    opt.loadfile (1,1) logical = true
    opt.overwrite (1,1) logical = false
end
persistent anovameta
if isempty(anovameta)
    if opt.loadfile
        load('anovameta.mat','anovameta')
    else
%         meta=ephys.util.load_meta();
        [~,~,sessmap]=ephys.sessid2path(0);
        homedir=ephys.util.getHomedir('type','raw');
        anovameta=struct();
        [anovameta.sess,anovameta.allcid,anovameta.anovap,anovameta.dur_selidx]=deal([]);
        for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
            disp(sesskey)
            fpath=fullfile(homedir,sessmap(sesskey),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');
            suid=h5read(fpath,'/SU_id');

            sensible_sel=ismember(trials(:,8),[3,6]) & all(trials(:,9:10)~=0,2) & ismember(trials(:,5),[4,8]);

            trials=trials(sensible_sel,:);
            fr=fr(sensible_sel,:,:);
%1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
            anovanps=nan(size(fr,2),7); 

            for su=1:size(fr,2)
                frmat=squeeze(fr(:,su,4:7));
                sufrvec=frmat(:);
                sampvec=repmat(trials(:,5),4,1);
                durvec=repmat(trials(:,8),4,1);
                binvec=reshape(repmat((1:4),size(fr,1),1),[],1);
                anovanps(su,:)=anovan(sufrvec,{sampvec,durvec,binvec},'model','full','display','off');

            end

            d=mean(sum(fr(trials(:,8)==3,:,4:7),3))-mean(sum(fr(trials(:,8)==6,:,4:7),3));
            s=mean(sum(fr(trials(:,8)==3,:,4:7),3))+mean(sum(fr(trials(:,8)==6,:,4:7),3));
             s(s==0)=-1;
            anovameta.dur_selidx=[anovameta.dur_selidx;(d./s).'];
            anovameta.sess=[anovameta.sess;repmat(sesskey,size(fr,2),1)];
            anovameta.allcid=[anovameta.allcid;suid];
            anovameta.anovap=[anovameta.anovap;anovanps];
        end
    end
end
anovameta_=anovameta;
if opt.overwrite
    save('anovameta.mat','anovameta')
end