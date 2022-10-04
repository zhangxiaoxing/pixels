function kwhmeta_=selectivity_kwh(opt)
arguments
    opt.loadfile (1,1) logical = false
    opt.overwrite (1,1) logical = false
    opt.per_bin (1,1) logical = true
end

persistent kwh_meta opt_
if isempty(kwh_meta) || ~isequaln(opt,opt_)
    if opt.loadfile
        load('kwh_meta.mat','kwh_meta')
    else
        [~,~,sessmap]=ephys.sessid2path(0);
        homedir=ephys.util.getHomedir('type','raw');
        kwh_meta=struct();
        [kwh_meta.sess,kwh_meta.allcid,kwh_meta.kwh_p,kwh_meta.fr,...
            kwh_meta.fdr_sel,kwh_meta.b,kwh_meta.in_m]=deal([]);
        for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
            disp(sesskey)
            fpath=fullfile(homedir,sessmap(sesskey),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');
            suid=h5read(fpath,'/SU_id');

            valid_sel=ismember(trials(:,8),[3,6]) & all(trials(:,9:10)~=0,2) & ismember(trials(:,5),[4,8]);

            trials=trials(valid_sel,:);
            fr=fr(valid_sel,:,:);

            if opt.per_bin
                %kruskal_wallis_p * 3 bins
                kwh_ps=nan(size(fr,2),3);
                fdr_sel=nan(size(fr,2),3);
                b=nan(size(fr,2),3,2);
                in_m=nan(size(fr,2),3,2);
            else
                kwh_ps=nan(size(fr,2),1);
            end
            if ~opt.per_bin
                for su=1:size(fr,2)
                    frmat=squeeze(mean(fr(:,su,5:7),3));
                    sufrvec=frmat;
                    sampvec=trials(:,5);
                    durvec=trials(:,8);
                    kwh_ps(su,:)=kruskalwallis(sufrvec,sampvec*10+durvec,'off');
                end
            else %ANOVA3
                for su=1:size(fr,2)
                    for bin=5:7
                        frmat=squeeze(fr(:,su,bin));
                        sufrvec=frmat(:);
                        sampvec=trials(:,5);
                        durvec=trials(:,8);
                        kwh_ps(su,bin-4)=kruskalwallis(sufrvec,sampvec*10+durvec,'off');
%                         [b(su,bin-4,:),~,~,in_m(su,bin-4,:)]=stepwisefit([sampvec==4,durvec==3],sufrvec,"display","off");
                        mdl=stepwiselm([sampvec,durvec],sufrvec,'interactions','ResponseVar','FR','PredictorVars',{'Sample','Duration'},...
                'CategoricalVars',{'Sample','Duration'},'Verbose',0);
                        if mdl.NumCoefficients>2
                            keyboard();
                        end
                    end
                    fdr_sel(su,:)=mafdr(squeeze(kwh_ps(su,1,:)),'BHFDR',true);
                end
            end

            s1d3=squeeze(mean(fr(trials(:,5)==4 & trials(:,8)==3,:,5:7)));
            s2d3=squeeze(mean(fr(trials(:,5)==8 & trials(:,8)==3,:,5:7)));
            s1d6=squeeze(mean(fr(trials(:,5)==4 & trials(:,8)==6,:,5:7)));
            s2d6=squeeze(mean(fr(trials(:,5)==8 & trials(:,8)==6,:,5:7)));

            kwh_meta.fr=[kwh_meta.fr;s1d3,s2d3,s1d6,s2d6];
            kwh_meta.sess=[kwh_meta.sess;repmat(sesskey,size(fr,2),1)];
            kwh_meta.allcid=[kwh_meta.allcid;suid];
            kwh_meta.kwh_p=[kwh_meta.kwh_p;kwh_ps];
            kwh_meta.b=[kwh_meta.b;b];
            kwh_meta.in_m=[kwh_meta.in_m;in_m];

            if exist('fdr_sel','var')
                kwh_meta.fdr_sel   =[kwh_meta.fdr_sel   ;fdr_sel   ];
            end
        end
        if ~opt.per_bin
            kwh_meta.selective    =kwh_meta.kwh_p(:,1)<0.05;
        else %per_bin
            kwh_meta.selective    =any(kwh_meta.fdr_sel<0.05,2);
        end
        kwh_meta.non_mem = ~kwh_meta.selective;
    end
end
opt_=opt;
kwhmeta_=kwh_meta;
if opt.overwrite
    save('kwh_meta.mat','kwh_meta')
end