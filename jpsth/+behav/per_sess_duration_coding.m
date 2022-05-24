function out=per_sess_duration_coding(opt)
arguments
    opt.new_data (1,1) logical=true
    opt.calc_dec (1,1) logical=true
    opt.plot_dec (1,1) logical=true
end
%% gen data

if opt.new_data

    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta();
    homedir=ephys.util.getHomedir('type','raw');
    for ii=54%reshape(cell2mat(sessmap.keys()),1,[])
        disp(ii)
        fpath=fullfile(homedir,sessmap(ii),"FR_All_1000.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
%         suid=h5read(fpath,'/SU_id');

        sesssel=meta.sess==ii;
        regsel=ismember(meta.reg_tree(1,sesssel),{'CH','BS'});
        if ~any(regsel)
            continue
        end
        fr=fr(:,regsel,:);

        fr_t_align=nan(size(fr,1),size(fr,2)); %trial, su, pre-test bin-averaged
%         fr_t_align(trials(:,8)==3,:)=mean(fr(trials(:,8)==3,:,23:28),3); % late delay
%         fr_t_align(trials(:,8)==6,:)=mean(fr(trials(:,8)==6,:,35:40),3); % late delay
        fr_t_align(trials(:,8)==3,:)=mean(fr(trials(:,8)==3,:,5:7),3); % late delay
        fr_t_align(trials(:,8)==6,:)=mean(fr(trials(:,8)==6,:,5:7),3); % late delay
        fr_t_align=normalize(fr_t_align,1,'range');
        c6sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==6 & ismember(trials(:,5),[4 8]);
        c3sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==3 & ismember(trials(:,5),[4 8]);
        e6sel=trials(:,10)==0 & trials(:,8)==6 & ismember(trials(:,5),[4 8]);
        e3sel=trials(:,10)==0 & trials(:,8)==3 & ismember(trials(:,5),[4 8]);

        cS1sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==4 & ismember(trials(:,8),[3 6]);
        cS2sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==8 & ismember(trials(:,8),[3 6]);
        eS1sel=trials(:,10)==0 & trials(:,5)==4 & ismember(trials(:,8),[3 6]);
        eS2sel=trials(:,10)==0 & trials(:,5)==8 & ismember(trials(:,8),[3 6]);
%%%%%%%%%%%%%%%%%%%%%%
        dur_vec=mean(fr_t_align(c6sel,:))-mean(fr_t_align(c3sel,:));
        dur_cd=dur_vec/norm(dur_vec);
        
        c6p=arrayfun(@(x) fr_t_align(x,:)*dur_cd.',find(c6sel));
        c3p=arrayfun(@(x) fr_t_align(x,:)*dur_cd.',find(c3sel));

        e6p=arrayfun(@(x) fr_t_align(x,:)*dur_cd.',find(e6sel));
        e3p=arrayfun(@(x) fr_t_align(x,:)*dur_cd.',find(e3sel));
%%%%%%%%%%%%%%%%%%
        olf_vec=mean(fr_t_align(cS1sel,:))-mean(fr_t_align(cS2sel,:));
        olf_cd=olf_vec/norm(olf_vec);
        
        cS1p=arrayfun(@(x) fr_t_align(x,:)*olf_cd.',find(cS1sel));
        cS2p=arrayfun(@(x) fr_t_align(x,:)*olf_cd.',find(cS2sel));

        eS1p=arrayfun(@(x) fr_t_align(x,:)*olf_cd.',find(eS1sel));
        eS2p=arrayfun(@(x) fr_t_align(x,:)*olf_cd.',find(eS2sel));
%%%%%%%%%%%%%%%%%%%

        fh=figure('Color','w','Position',[32,32,500,225]);
        subplot(1,2,2)
        hold on;
        swarmchart(ones(size(c6p)),c6p,  4,'r','filled','o','MarkerFaceAlpha','0.2','MarkerEdgeColor','none')
        swarmchart(2*ones(size(c3p)),c3p,4,'b','filled','o','MarkerFaceAlpha','0.2','MarkerEdgeColor','none')
        swarmchart(3*ones(size(e6p)),e6p,4,'r','filled','o','MarkerFaceAlpha','0.2','MarkerEdgeColor','none')
        swarmchart(4*ones(size(e3p)),e3p,4,'b','filled','o','MarkerFaceAlpha','0.2','MarkerEdgeColor','none')

        boxchart(ones(size(c6p)),c6p,'LineWidth',1,'BoxFaceColor','k','WhiskerLineColor','k','BoxFaceAlpha',0,'MarkerStyle','none')
        boxchart(2*ones(size(c3p)),c3p,'LineWidth',1,'BoxFaceColor','k','WhiskerLineColor','k','BoxFaceAlpha',0,'MarkerStyle','none')
        boxchart(3*ones(size(e6p)),e6p,'LineWidth',1,'BoxFaceColor','k','WhiskerLineColor','k','BoxFaceAlpha',0,'MarkerStyle','none')
        boxchart(4*ones(size(e3p)),e3p,'LineWidth',1,'BoxFaceColor','k','WhiskerLineColor','k','BoxFaceAlpha',0,'MarkerStyle','none')

        ylabel('Duration projection (a.u.)');
        set(gca(),'XTick',1:4,'XTickLabel',{'6s','3s','6s','3s'});

        subplot(1,2,1)
        hold on;
        swarmchart(ones(size(cS1p)),cS1p,  4,'r','filled','o','MarkerFaceAlpha','0.2','MarkerEdgeColor','none')
        swarmchart(2*ones(size(cS2p)),cS2p,4,'b','filled','o','MarkerFaceAlpha','0.2','MarkerEdgeColor','none')
        swarmchart(3*ones(size(eS1p)),eS1p,4,'r','filled','o','MarkerFaceAlpha','0.2','MarkerEdgeColor','none')
        swarmchart(4*ones(size(eS2p)),eS2p,4,'b','filled','o','MarkerFaceAlpha','0.2','MarkerEdgeColor','none')
        
        boxchart(ones(size(cS1p)),cS1p,'LineWidth',1,'BoxFaceColor','k','WhiskerLineColor','k','BoxFaceAlpha',0,'MarkerStyle','none')
        boxchart(2*ones(size(cS2p)),cS2p,'LineWidth',1,'BoxFaceColor','k','WhiskerLineColor','k','BoxFaceAlpha',0,'MarkerStyle','none')
        boxchart(3*ones(size(eS1p)),eS1p,'LineWidth',1,'BoxFaceColor','k','WhiskerLineColor','k','BoxFaceAlpha',0,'MarkerStyle','none')
        boxchart(4*ones(size(eS2p)),eS2p,'LineWidth',1,'BoxFaceColor','k','WhiskerLineColor','k','BoxFaceAlpha',0,'MarkerStyle','none')

        ylabel('Sensory projection (a.u.)');
        set(gca(),'XTick',1:4,'XTickLabel',{'S1','S2','S1','S2'});
        sgtitle(ii);
        exportgraphics(fh,'sens_dur_cd_proj_sess_sc.pdf','ContentType','vector');
        waitfor(fh);
    end
end
