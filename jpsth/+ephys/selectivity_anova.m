function anovameta_=selectivity_anova(opt)
arguments
    opt.loadfile (1,1) logical = false
    opt.overwrite (1,1) logical = false
    opt.merge_time_bin (1,1) logical = false
end
persistent anovameta opt_
if isempty(anovameta) || ~isequaln(opt,opt_)
    if opt.loadfile
        load('anovameta.mat','anovameta')
    else
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

            if opt.merge_time_bin
                %1:Samp,2:Dur,3:Samp*Dur
                anovanps=nan(size(fr,2),3);
            else
                %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
                anovanps=nan(size(fr,2),7);
            end
            if opt.merge_time_bin
                for su=1:size(fr,2)
                    frmat=squeeze(mean(fr(:,su,5:7),3));
                    sufrvec=frmat;
                    sampvec=trials(:,5);
                    durvec=trials(:,8);
                    anovanps(su,:)=anovan(sufrvec,{sampvec,durvec},'model','full','display','off');
                end
            else
                for su=1:size(fr,2)
                    frmat=squeeze(fr(:,su,5:7));
                    sufrvec=frmat(:);
                    sampvec=repmat(trials(:,5),3,1);
                    durvec=repmat(trials(:,8),3,1);
                    binvec=reshape(repmat((1:3),size(fr,1),1),[],1);
                    anovanps(su,:)=anovan(sufrvec,{sampvec,durvec,binvec},'model','full','display','off');
                end
            end

            s1d3=mean(sum(fr(trials(:,5)==4 & trials(:,8)==3,:,4:7),3));
            s2d3=mean(sum(fr(trials(:,5)==8 & trials(:,8)==3,:,4:7),3));
            s1d6=mean(sum(fr(trials(:,5)==4 & trials(:,8)==6,:,4:7),3));
            s2d6=mean(sum(fr(trials(:,5)==8 & trials(:,8)==6,:,4:7),3));

            anovameta.dur_selidx=[anovameta.dur_selidx;s1d3.',s2d3.',s1d6.',s2d6.'];
            anovameta.sess=[anovameta.sess;repmat(sesskey,size(fr,2),1)];
            anovameta.allcid=[anovameta.allcid;suid];
            anovameta.anovap=[anovameta.anovap;anovanps];
        end
        if opt.merge_time_bin
            anovameta.sens_only=anovameta.anovap(:,1)<0.05 & ~any(anovameta.anovap(:,2:3)<0.05,2);
            anovameta.dur_only=anovameta.anovap(:,2)<0.05 & ~any(anovameta.anovap(:,[1 3])<0.05,2);
            anovameta.mixed=all(anovameta.anovap(:,1:2)<0.05,2) | anovameta.anovap(:,3)<0.05;
            anovameta.non_mem= true(size(anovameta.allcid)) & ~(anovameta.sens_only | anovameta.dur_only | anovameta.mixed);
        end
    end
end
opt_=opt;
anovameta_=anovameta;
if opt.overwrite
    save('anovameta.mat','anovameta')
end