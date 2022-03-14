%%CONST
meta=ephys.util.load_meta();
[~,~,sessmap]=ephys.sessid2path(0);
homedir=ephys.util.getHomedir('type','raw');

epochs=["ITI","sample","early_delay","late_delay","test"];
%%
map_cells=cell(5,6);
for epochi=1:5
    epoch=epochs(epochi);
    anovameta=struct();
    [anovameta.sess,anovameta.allcid,anovameta.anovap]=deal([]);

    for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
        disp(sesskey)
        fpath=fullfile(homedir,sessmap(sesskey),"FR_All_ 250.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        dur_resp=behav.tag_block(trials,'wt',false);
        block_meta=[trials,dur_resp(:,end)];

        sensible_sel=(trials(:,8)==3 | trials(:,8)==6) & all(trials(:,9:10)~=0,2);
        trials=trials(sensible_sel,:);
        fr=fr(sensible_sel,:,:);
        block_meta=block_meta(sensible_sel,:);

        anovanps=nan(size(fr,2),7); %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin

        for su=1:size(fr,2)
            switch epoch
                case 'late_delay'
                    frmat=squeeze(fr(:,su,23:28));
                    frmat(trials(:,8)==6,:)=squeeze(fr(trials(:,8)==6,su,35:40));
                    frvec=sum(frmat./4,2)./1.5;
                case 'sample'
                    frmat=squeeze(fr(:,su,13:16));
                    frvec=sum(frmat./4,2);
                case 'ITI'
                    frmat=squeeze(fr(:,su,3:10));
                    frvec=sum(frmat./4,2)./2;
                case 'early_delay'
                    frmat=squeeze(fr(:,su,17:22));
                    frvec=sum(frmat./4,2)./1.5;
                case 'test'
                    frmat=squeeze(fr(:,su,29:32));
                    frmat(trials(:,8)==6,:)=squeeze(fr(trials(:,8)==6,su,41:44));
                    frvec=sum(frmat./4,2)./1.5;
            end
            sampvec=trials(:,5);
            durvec=trials(:,8);
            seqvec=block_meta(:,11);
            if strcmp(epoch,'ITI')
                anovanps(su,[2 3 6])=anovan(frvec,{durvec,seqvec},'model','full','display','off');
            else
                anovanps(su,:)=anovan(frvec,{sampvec,durvec,seqvec},'model','full','display','off');
            end
        end

        anovameta.sess=[anovameta.sess;repmat(sesskey,size(fr,2),1)];
        anovameta.allcid=[anovameta.allcid;suid];
        anovameta.anovap=[anovameta.anovap;anovanps];
    end


    %% scatter common
    meta=ephys.util.load_meta();
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));

    BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
    CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
    grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));

    cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
    grey_regs=grey_regs(cnt>100);

    %% samp v dur v seq

    samp_any_sel=any(anovameta.anovap(:,[1,4,5,7])<0.05,2);
    dur_any_sel=any(anovameta.anovap(:,[2,4,6,7])<0.05,2);
    seq_any_sel=any(anovameta.anovap(:,[3,5,6,7])<0.05,2);

    samp_only_sel=samp_any_sel & ~dur_any_sel & ~seq_any_sel;
    dur_only_sel=dur_any_sel & ~samp_any_sel & ~seq_any_sel;
    seq_only_sel=seq_any_sel & ~samp_any_sel & ~dur_any_sel;

    samp_any_reg_map=containers.Map();
    samp_only_reg_map=containers.Map();

    dur_any_reg_map=containers.Map();
    dur_only_reg_map=containers.Map();

    seq_any_reg_map=containers.Map();
    seq_only_reg_map=containers.Map();

    for r=grey_regs
        cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
        samp_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & samp_any_sel);
        samp_only_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & samp_only_sel);

        dur_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & dur_any_sel);
        dur_only_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & dur_only_sel);

        seq_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & seq_any_sel);
        seq_only_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & seq_only_sel);

        samp_any_reg_map(r{1})=[samp_any_cnt/cnt,samp_any_cnt,cnt];
        samp_only_reg_map(r{1})=[samp_only_cnt/cnt,samp_only_cnt,cnt];

        dur_any_reg_map(r{1})=[dur_any_cnt/cnt,dur_any_cnt,cnt];
        dur_only_reg_map(r{1})=[dur_only_cnt/cnt,dur_only_cnt,cnt];

        seq_any_reg_map(r{1})=[seq_any_cnt/cnt,seq_any_cnt,cnt];
        seq_only_reg_map(r{1})=[seq_only_cnt/cnt,seq_only_cnt,cnt];
    end

    map_cells(epochi,:)={samp_any_reg_map,samp_only_reg_map,...
        dur_any_reg_map,dur_only_reg_map,...
        seq_any_reg_map,seq_only_reg_map};

    maps={samp_any_reg_map,dur_any_reg_map;...
        samp_any_reg_map,seq_any_reg_map;...
        dur_any_reg_map,seq_any_reg_map;...
        samp_only_reg_map,dur_only_reg_map;...
        samp_only_reg_map,seq_only_reg_map;...
        dur_only_reg_map,seq_only_reg_map...
        };

    lbls={'Mixed sample coding proportion','Mixed duration coding proportion';...
        'Mixed sample coding proportion','Mixed block-position coding proportion';...
        'Mixed duration coding proportion','Mixed block-position coding proportion';...
        'Exclusive sample coding proportion','Exclusive duration coding proportion';...
        'Exclusive sample coding proportion','Exclusive block-position coding proportion';...
        'Exclusive duration coding proportion','Exclusive block-position coding proportion'};

    figure('Color','w','Position',[32,32,1200,800])
    for ii=1:6
        if strcmp(epoch,'ITI') && ~ismember(ii,[3,6])
            continue
        end
        xmap=maps{ii,1};%
        ymap=maps{ii,2};%

        x_mat=cell2mat(xmap.values(grey_regs).');
        y_mat=cell2mat(ymap.values(grey_regs).');

        [r,p]=corr(x_mat(:,1),y_mat(:,1));

        subplot(2,3,ii);
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


%% cross-time correlation heatmap
% any selective, [1 3 5]; exclusive selective, [2 4 6]
% samp:samp, dur:dur, seq:seq
for stype=["Any","Exclusive"]
rmat=cell(2,3);
cols=1:2:5;
if strcmp(stype,'Exclusive')
    cols=cols+1;
end
for colidx=1:3
    rmat{1,colidx}=nan(5,5);
    for ii=1:5
        for jj=ii:5
            if colidx==1 && (ii==1 || jj==1)
                continue
            elseif ii==jj
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

% samp:dur, samp:seq, dur:seq
col_vs=[1,3;1,5;3,5];
if strcmp(stype,'Exclusive')
    col_vs=col_vs+1;
end
for vs_idx=1:3
    rmat{2,vs_idx}=nan(5,5);
    for ii=1:5
        for jj=1:5
            if col_vs(vs_idx,1)<3 && ii==1
                continue
            else
                x_mat=cell2mat(map_cells{ii,col_vs(vs_idx,1)}.values(grey_regs).');
                y_mat=cell2mat(map_cells{jj,col_vs(vs_idx,2)}.values(grey_regs).');
                r=corr(x_mat(:,1),y_mat(:,1));
                rmat{2,vs_idx}(ii,jj)=r;
            end
        end
    end
end

axlbl={'Sample','Sample';...
    'Duration','Duration';...
    'Block-position','Block-position';...
    'Sample','Duration';...
    'Sample','Block-position';...
    'Duration','Block-position'};

figure('Color','w','Position',[32,32,1200,800])
for ii=1:6
subplot(2,3,ii)
imagesc(rmat{(ii>3)+1,ii-(ii>3)*3},[-1,1]);
colormap('turbo')
ticklbl={'ITI','SMP','EDL','LDL','TST'};
set(gca,'YDir','normal','XTick',1:5,'XTickLabel',ticklbl,'YTick',1:5,'YTickLabel',ticklbl)
xlabel(axlbl{ii,1});
ylabel(axlbl{ii,2});
end
sgtitle(stype)
end

%% cross-time correlation scatter
% any selective, [1 3 5]; exclusive selective, [2 4 6]
% samp:samp, dur:dur, seq:seq
figure('Color','w','Position',[32,32,1200,800])
for stype=["Any","Exclusive"]
    rmat=cell(2,3);
    cols=1:2:5;
    if strcmp(stype,'Exclusive')
        cols=cols+1;
    end
    lbls={'Sample early delay', 'Sample late delay';...
        'Duration early delay', 'Duration late delay';...
        'Block-position early delay', 'Block-position late delay'};

    for colidx=1:3
        x_mat=cell2mat(map_cells{3,cols(colidx)}.values(grey_regs).');
        y_mat=cell2mat(map_cells{4,cols(colidx)}.values(grey_regs).');
        r=corr(x_mat(:,1),y_mat(:,1));
        if strcmp(stype,'Any')
            subplot(2,3,colidx);
        else
            subplot(2,3,colidx+3)
        end
        hold on
        for rr=reshape(grey_regs,1,[])
            c=ephys.getRegColor(rr{1},'large_area',true);
            xx=map_cells{3,cols(colidx)}(rr{1});
            yy=map_cells{4,cols(colidx)}(rr{1});
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


return

%%
