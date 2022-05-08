function anovameta_=selectivity_anova(opt)
arguments
    opt.loadfile (1,1) logical = false
    opt.overwrite (1,1) logical = false
    opt.merge_time_bin (1,1) logical = false
    opt.largest_varied_bin (1,1) logical = false
    opt.per_bin (1,1) logical = false
    opt.anova_model (1,:) char {mustBeMember(opt.anova_model,{'linear','full','handpick'})}='full'
end
assert(nnz([opt.merge_time_bin, opt.largest_varied_bin,opt.per_bin])<2,'Must choose one method')

if strcmp(opt.anova_model,'handpick')
    opt.anova_model=[1,0,0;0,1,0;0,0,1;1,1,0];
end

persistent anovameta opt_
if isempty(anovameta) || ~isequaln(opt,opt_)
    if opt.loadfile
        load('anovameta.mat','anovameta')
    else
        [~,~,sessmap]=ephys.sessid2path(0);
        homedir=ephys.util.getHomedir('type','raw');
        anovameta=struct();
        [anovameta.sess,anovameta.allcid,anovameta.anovap,anovameta.dur_selidx,...
            anovameta.fdr_sense,anovameta.fdr_dur,anovameta.fdr_interact]=deal([]);
        for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
            disp(sesskey)
            fpath=fullfile(homedir,sessmap(sesskey),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');
            suid=h5read(fpath,'/SU_id');

            sensible_sel=ismember(trials(:,8),[3,6]) & all(trials(:,9:10)~=0,2) & ismember(trials(:,5),[4,8]);

            trials=trials(sensible_sel,:);
            fr=fr(sensible_sel,:,:);

            if opt.merge_time_bin || opt.largest_varied_bin
                %1:Samp,2:Dur,3:Samp*Dur
                Ncoeff=2+strcmp(opt.anova_model,'full');
                anovanps=nan(size(fr,2),Ncoeff);
            elseif opt.per_bin
                %[1:Samp,2:Dur,3:Samp*Dur] * 3 bins
                anovanps=nan(size(fr,2),3,3);
                [fdr_sense,fdr_dur,fdr_interact]=deal(nan(size(fr,2),3));
            else
                if isa(opt.anova_model,'double')
                    anovanps=nan(size(fr,2),4);
                elseif strcmp(opt.anova_model,'full')
                    %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
                    anovanps=nan(size(fr,2),7);
                else
                    anovanps=nan(size(fr,2),3);
                end
            end
            if opt.merge_time_bin
                for su=1:size(fr,2)
                    frmat=squeeze(mean(fr(:,su,5:7),3));
                    sufrvec=frmat;
                    sampvec=trials(:,5);
                    durvec=trials(:,8);
                    anovanps(su,:)=anovan(sufrvec,{sampvec,durvec},'display','off','model',opt.anova_model);

                end
            elseif opt.largest_varied_bin
                for su=1:size(fr,2)
                    stds=std(squeeze(fr(:,su,5:7)));
                    [~,bin]=max(stds);
                    sufrvec=fr(:,su,bin+4);
                    sampvec=trials(:,5);
                    durvec=trials(:,8);
                    anovanps(su,:)=anovan(sufrvec,{sampvec,durvec},'display','off','model',opt.anova_model);
                end                
            elseif opt.per_bin %ANOVA3
                for su=1:size(fr,2)
                    for bin=5:7
                        frmat=squeeze(fr(:,su,bin));
                        sufrvec=frmat(:);
                        sampvec=trials(:,5);
                        durvec=trials(:,8);
                        anovanps(su,:,bin-4)=anovan(sufrvec,{sampvec,durvec},'model','full','display','off');
                    end
                    fdr_sense(su,:)=mafdr(squeeze(anovanps(su,1,:)),'BHFDR',true);
                    fdr_dur(su,:)=mafdr(squeeze(anovanps(su,2,:)),'BHFDR',true);
                    fdr_interact(su,:)=mafdr(squeeze(anovanps(su,3,:)),'BHFDR',true);
                end
            else
                for su=1:size(fr,2)
                    frmat=squeeze(fr(:,su,5:7));
                    sufrvec=frmat(:);
                    sampvec=repmat(trials(:,5),3,1);
                    durvec=repmat(trials(:,8),3,1);
                    binvec=reshape(repmat((1:3),size(fr,1),1),[],1);
                    
                    anovanps(su,:)=anovan(sufrvec,{sampvec,durvec,binvec},'model',opt.anova_model,'display','off');
                end
            end

            s1d3=squeeze(mean(fr(trials(:,5)==4 & trials(:,8)==3,:,5:7)));
            s2d3=squeeze(mean(fr(trials(:,5)==8 & trials(:,8)==3,:,5:7)));
            s1d6=squeeze(mean(fr(trials(:,5)==4 & trials(:,8)==6,:,5:7)));
            s2d6=squeeze(mean(fr(trials(:,5)==8 & trials(:,8)==6,:,5:7)));

            anovameta.dur_selidx=[anovameta.dur_selidx;s1d3,s2d3,s1d6,s2d6];
            anovameta.sess=[anovameta.sess;repmat(sesskey,size(fr,2),1)];
            anovameta.allcid=[anovameta.allcid;suid];
            anovameta.anovap=[anovameta.anovap;anovanps];
            if exist('fdr_sense','var')
                anovameta.fdr_sense   =[anovameta.fdr_sense   ;fdr_sense   ];
                anovameta.fdr_dur     =[anovameta.fdr_dur     ;fdr_dur     ];
                anovameta.fdr_interact=[anovameta.fdr_interact;fdr_interact];
            end
        end
        if (opt.merge_time_bin || opt.largest_varied_bin)
            if strcmp(opt.anova_model,'full')
                anovameta.sens    =anovameta.anovap(:,1)<0.05;
                anovameta.dur     =anovameta.anovap(:,2)<0.05;
                anovameta.interact=anovameta.anovap(:,3)<0.05;
                anovameta.non_mem =true(size(anovameta.allcid)) & ~(anovameta.sens | anovameta.dur | anovameta.interact);
            else
                anovameta.sens    =anovameta.anovap(:,1)<0.05;
                anovameta.dur     =anovameta.anovap(:,2)<0.05;
                anovameta.interact=false(size(anovameta.allcid));
                anovameta.non_mem =true(size(anovameta.allcid)) & ~(anovameta.sens | anovameta.dur | anovameta.interact);
             end
        elseif opt.per_bin
                anovameta.sens    =any(anovameta.fdr_sense<0.05,2);
                anovameta.dur     =any(anovameta.fdr_dur<0.05,2);
                anovameta.interact=any(anovameta.fdr_interact<0.05,2);
                anovameta.non_mem =true(size(anovameta.allcid)) & ~(anovameta.sens | anovameta.dur | anovameta.interact);
        else
            if ~strcmp(opt.anova_model,'linear') 
            %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
                anovameta.sens    =anovameta.anovap(:,1)<0.05;
                anovameta.dur     =anovameta.anovap(:,2)<0.05;
                anovameta.time_bin=anovameta.anovap(:,3)<0.05;
                anovameta.interact=anovameta.anovap(:,4)<0.05;
                anovameta.non_mem =true(size(anovameta.allcid)) & ~(anovameta.sens | anovameta.dur | anovameta.interact);
            else
                anovameta.sens    =anovameta.anovap(:,1)<0.05;
                anovameta.dur     =anovameta.anovap(:,2)<0.05;
                anovameta.time_bin=anovameta.anovap(:,3)<0.05;
                anovameta.interact=false(size(anovameta.allcid));
                anovameta.non_mem =true(size(anovameta.allcid)) & ~(anovameta.sens | anovameta.dur | anovameta.interact);
            end
        end
    end
end
opt_=opt;
anovameta_=anovameta;
if opt.overwrite
    save('anovameta.mat','anovameta')
end