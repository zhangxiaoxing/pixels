%%CONST
meta=ephys.util.load_meta();
[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir('type','raw');

% epochs=["sample","early_delay","late_delay","test","postreward","presample"];
epochs=["alt_delay"];
%%
if ~exist('cross_ep_anovameta','var')
    cross_ep_anovameta=struct();
end
map_cells=cell(numel(epochs),4);
for epochi=1:(numel(epochs))
    epoch=epochs(epochi);
    anovameta=struct();
    [anovameta.sess,anovameta.allcid,anovameta.anovap]=deal([]);

    for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
        disp(sesskey)
        fpath=fullfile(homedir,sessmap(sesskey),"FR_All_ 250.hdf5");
        if strcmp(epoch,'postreward') || strcmp(epoch,'presample')
            [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sesskey,"keep_trial",true,"align_test",true);
            cfg=struct();
            cfg.binsize=0.25;
            cfg.keeptrials='yes';
            FT_PSTH=ft_spike_psth(cfg, FT_SPIKE);
            fr=FT_PSTH.trial;
        else
            fr=h5read(fpath,'/FR_All');
        end

        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        sensible_sel=ismember(trials(:,8),[3,6]) & ismember(trials(:,5),[4,8]) & all(trials(:,9:10)~=0,2);
        trials=trials(sensible_sel,:);
        fr=fr(sensible_sel,:,:);

        anovanps=nan(size(fr,2),3); %1:Samp,2:Dur,3:Samp*Dur

        for su=1:size(fr,2)
            switch epoch
                case 'late_delay'
                    frmat=squeeze(fr(:,su,23:28));
                    frmat(trials(:,8)==6,:)=squeeze(fr(trials(:,8)==6,su,35:40));
                    frvec=sum(frmat./4,2)./1.5;
                case 'sample'
                    frmat=squeeze(fr(:,su,13:16));
                    frvec=sum(frmat./4,2);
                case 'presample'
                    frmat=squeeze(fr(:,su,37:56));
                    frvec=sum(frmat./4,2)./5;
                case 'early_delay'
                    frmat=squeeze(fr(:,su,17:22));
                    frvec=sum(frmat./4,2)./1.5;
                case 'test'
                    frmat=squeeze(fr(:,su,29:32));
                    frmat(trials(:,8)==6,:)=squeeze(fr(trials(:,8)==6,su,41:44));
                    frvec=sum(frmat./4,2);
                case 'postreward'
                    frmat=squeeze(fr(:,su,17:36));
                    frvec=sum(frmat./4,2)./5;

                case 'alt_delay'
                    frmat=squeeze(fr(:,su,17:28));
                    frvec=sum(frmat./4,2)./3;
            end
            sampvec=trials(:,5);
            durvec=trials(:,8);
            if mean(frvec<0.5)
                anovanps(su,:)=ones(1,3);
            else
                anovanps(su,:)=anovan(frvec,{sampvec,durvec},'model','full','continuous',[],'display','off');
            end
        end

        anovameta.sess=[anovameta.sess;repmat(sesskey,size(fr,2),1)];
        anovameta.allcid=[anovameta.allcid;suid];
        anovameta.anovap=[anovameta.anovap;anovanps];
    end
    
    if false % plot
        %% scatter common
        idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
        grey_regs=ephys.getGreyRegs();
        %% samp v dur v seq

 

        map_cells(epochi,:)={samp_any_reg_map,samp_only_reg_map,...
            dur_any_reg_map,dur_only_reg_map};

        %     if false % read back from previous save
        %         samp_any_reg_map=map_cells{1,1};
        %         samp_only_reg_map=map_cells{1,2};
        %         dur_any_reg_map=map_cells{1,3};
        %         dur_only_reg_map=map_cells{1,4};
        %         seq_any_reg_map=map_cells{1,5};
        %         seq_only_reg_map=map_cells{1,6};
        %     end

        maps={samp_any_reg_map,dur_any_reg_map;...
            samp_only_reg_map,dur_only_reg_map};

        lbls={'Mixed sample coding proportion','Mixed duration coding proportion';...
            'Exclusive sample coding proportion','Exclusive duration coding proportion'};

        figure('Color','w','Position',[32,32,900,400])
        for ii=1:2
            xmap=maps{ii,1};%
            ymap=maps{ii,2};%

            x_mat=cell2mat(xmap.values(grey_regs).');
            y_mat=cell2mat(ymap.values(grey_regs).');

            [r,p]=corr(x_mat(:,1),y_mat(:,1));

            subplot(1,2,ii);
            hold on
            for rr=reshape(grey_regs,1,[])
                c=ephys.getRegColor(rr{1},'large_area',true);
                xx=xmap(rr{1});
                yy=ymap(rr{1});
                scatter(xx(1),yy(1),9,c,'filled','o')
                text(xx(1),yy(1),rr{1},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
            end
            title(sprintf('r=%.2f,p=%.2f',r,p),'FontSize',10);
            xlabel(lbls{ii,1});
            ylabel(lbls{ii,2});
        end
        sgtitle(epoch)
    end
    cross_ep_anovameta.(epoch)=anovameta;
end
save('cross_ep_anovameta.mat','cross_ep_anovameta');
%% cross-time correlation heatmap
% any selective, [1 3]; exclusive selective, [2 4]
% samp:samp, dur:dur
for stype=["Mixed","Exclusive"]
    rmat=cell(1,3);
    cols=1:2:3;
    if strcmp(stype,'Exclusive')
        cols=cols+1;
    end
    for colidx=1:2
        rmat{1,colidx}=nan(6,6);
        for ii=1:6
            for jj=ii:6
                %             if colidx==1 && (ii==1 || jj==1)
                %                 continue
                %             else
                if ii==jj
                    [rmat{1,colidx}(ii,jj),rmat{1,colidx}(jj,ii)]=deal(1);
                else
                    x_mat=cell2mat(map_cells{ii,cols(colidx)}.values(grey_regs).');
                    y_mat=cell2mat(map_cells{jj,cols(colidx)}.values(grey_regs).');
                    r=corr(x_mat(:,1),y_mat(:,1));
                    [rmat{1,colidx}(ii,jj),rmat{1,colidx}(jj,ii)]=deal(r);
                end
            end
        end
    end

    % samp:dur
    col_vs=[1,3];
    if strcmp(stype,'Exclusive')
        col_vs=col_vs+1;
    end

    rmat{1,3}=nan(6,6);
    for ii=1:6
        for jj=1:6
            x_mat=cell2mat(map_cells{ii,col_vs(1)}.values(grey_regs).');
            y_mat=cell2mat(map_cells{jj,col_vs(2)}.values(grey_regs).');
            r=corr(x_mat(:,1),y_mat(:,1));
            rmat{1,3}(ii,jj)=r;
        end
    end


    axlbl={'Sample','Sample';...
        'Duration','Duration';...
        'Sample','Duration'};

    figure('Color','w','Position',[32,32,1024,300])
    for ii=1:3
        subplot(1,3,ii)
        imagesc(rmat{1,ii},[-1,1]);
        colormap('turbo')
        ticklbl={'SMP','EDL','LDL','TST','PosT','PreS'};
        set(gca,'YDir','normal','XTick',1:6,'XTickLabel',ticklbl,'YTick',1:6,'YTickLabel',ticklbl)
        xlabel(axlbl{ii,1});
        ylabel(axlbl{ii,2});
    end
    sgtitle(stype)

end

%% cross-time correlation scatter
% any selective, [1 3 5]; exclusive selective, [2 4 6]
% samp:samp, dur:dur, seq:seq
figure('Color','w','Position',[32,32,800,400])
for stype=["Mixed","Exclusive"]
    rmat=cell(2,3);
    cols=1:2:5;
    if strcmp(stype,'Exclusive')
        cols=cols+1;
    end
    lbls={'Sample early delay', 'Sample late delay';...
        'Duration early delay', 'Duration late delay';...
        'Block-position early delay', 'Block-position late delay'};

    for colidx=1:3
        x_mat=cell2mat(map_cells{2,cols(colidx)}.values(grey_regs).');
        y_mat=cell2mat(map_cells{3,cols(colidx)}.values(grey_regs).');
        r=corr(x_mat(:,1),y_mat(:,1));
        if strcmp(stype,'Mixed')
            subplot(2,3,colidx);
        else
            subplot(2,3,colidx+3)
        end
        hold on
        for rr=reshape(grey_regs,1,[])
            c=ephys.getRegColor(rr{1},'large_area',true);
            xx=map_cells{2,cols(colidx)}(rr{1});
            yy=map_cells{3,cols(colidx)}(rr{1});
            scatter(xx(1),yy(1),9,c,'filled','o')
            text(xx(1),yy(1),rr{1},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
        end
        title(sprintf('r=%.2f,p=%.2f',r,p),'FontSize',10);
        maxspan=max([xlim(),ylim()]);
        xlim([0,maxspan]);
        ylim([0,maxspan]);
        plot([0,maxspan],[0,maxspan],'--k')

        xlabel(lbls{colidx,1});
        ylabel(lbls{colidx,2});
    end
end
sgtitle('Early- vs late-delay, mixed selective (top) or exclusively selective (bottom)')

figure('Color','w','Position',[32,32,600,300])
pathstr=strjoin({mfilename('fullpath'),matlab.desktop.editor.getActiveFilename},'|');
text(0.5,0.5,pathstr,'HorizontalAlignment','center','VerticalAlignment','middle')


%% map_cells (epoch x sample-dur-pos) -> GLM


return

