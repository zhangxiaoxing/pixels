%TODO brain region filter, olfaction filter.
function [fh,out]=pct_decoding(p_list,bound,opt)
arguments
    p_list
    bound
    opt.new_data (1,1) logical=true
    %     opt.plot_PCA (1,1) logical=false
    opt.calc_dec (1,1) logical=true
    opt.plot_dec (1,1) logical=true
    opt.ranksum_stats (1,1) logical =false
    opt.lblidx (1,1) double {mustBeMember(opt.lblidx,[5,8])} = 8 %5 for sample 8 for duration
    opt.rpt (1,1) double {mustBeInteger,mustBePositive} = 100
    opt.n_su (1,:) double {mustBeInteger,mustBePositive} =12:12:48;
    opt.cmap (1,:) char {mustBeMember(opt.cmap,{'parula','cool'})}
    opt.cross (1,1) logical = false
end

if opt.lblidx==5
    behav_tags=[4,8];
    dec_tag="olf";
else
    behav_tags=[3,6];
    dec_tag="dur";
end

if opt.cross
    dec_tag=dec_tag+"x";
end

ebound=reshape(bound,1,[]);
%% gen data

if opt.new_data
    decode_mat=struct();
    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
    homedir=ephys.util.getHomedir('type','raw');
    regsel=ismember(meta.reg_tree(1,:),{'CH','BS'}).';
    %% pct-loop
    frmap=struct();
    for bidx=1:(numel(ebound)-1)
        susel=p_list>ebound(bidx) & p_list<=ebound(bidx+1);
        for ii=reshape(cell2mat(sessmap.keys()),1,[])
            disp(ii)
            fpath=fullfile(homedir,sessmap(ii),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');
            sesssel=meta.sess==ii;
            fr_t_align=mean(fr(:,:,5:7),3);

            decode_mat.(sprintf('s%d',ii)).trials=trials;
            decode_mat.(sprintf('s%d',ii)).("band"+num2str(bidx))=fr_t_align(:,regsel(sesssel) & susel(sesssel));
        end
        frmap.(sprintf('fr%d',behav_tags(1))).("band"+num2str(bidx))=containers.Map('KeyType','char','ValueType','any');
        frmap.(sprintf('fr%d',behav_tags(2))).("band"+num2str(bidx))=containers.Map('KeyType','char','ValueType','any');
    end

    trlCount=[];
    for skey=reshape(fieldnames(decode_mat),1,[])
        trls=decode_mat.(skey{1}).trials;
        csel=trls(:,9)~=0 & trls(:,10)~=0 & ismember(trls(:,8),[3 6]) & ismember(trls(:,5),[4 8]);
        dec_c_lbl=trls(csel,opt.lblidx);
        %% TODO: more pct
        c_trlSess=[0,0];
        for behav_idx=1:2
            c_trlSess(behav_idx)=nnz(dec_c_lbl==behav_tags(behav_idx));
            for bidx=1:(numel(ebound)-1)
                cfrs=decode_mat.(skey{1}).("band"+num2str(bidx))(csel,:); % [nTrl,nSU]
                for suidx=1:size(cfrs,2) % TODO: vectorized batch init
                    frmap.(sprintf('fr%d',behav_tags(behav_idx))).("band"+num2str(bidx))(sprintf('%su%d',skey{1},suidx))=cfrs(dec_c_lbl==behav_tags(behav_idx),suidx);
                end
            end
        end
        trlCount=[trlCount;c_trlSess];
    end
    save(sprintf('pct_decode_fr_%s.mat',dec_tag),'frmap','trlCount');
else
    load(sprintf('pct_decode_fr_%s.mat',dec_tag),'frmap','trlCount');
end
%% available trials for least repeated condition
%figure();histogram(trlCount(:),1:32)


%% decoding
if opt.calc_dec
    min_trl=20;
    lbls=["fr"+num2str(behav_tags(1)),"fr"+num2str(behav_tags(2))];
    trl_cnt={trlCount};
    out=struct();

    for bidx=1:(numel(ebound)-1)
        disp("band #"+num2str(bidx))
        curr_keys=frmap.(char(lbls(1))).("band"+num2str(bidx)).keys();
        dur_trlSel=ismember(cellfun(@(x) str2double(regexp(x,'(?<=s)\d*(?=u)','match','once')),curr_keys),find(min(trl_cnt{1},[],2)>min_trl));
        curr_keys=curr_keys(dur_trlSel);

        for n_su=opt.n_su
            result=[];
            for resamp_rpt=1:opt.rpt%15
                sukeys=datasample(curr_keys,n_su,'replace',false);
                rawmat=[...
                    cellfun(@(x) min(x),frmap.(lbls(1)).("band"+num2str(bidx)).values(sukeys));...
                    cellfun(@(x) max(x),frmap.(lbls(1)).("band"+num2str(bidx)).values(sukeys));...
                    cellfun(@(x) min(x),frmap.(lbls(2)).("band"+num2str(bidx)).values(sukeys));...
                    cellfun(@(x) max(x),frmap.(lbls(2)).("band"+num2str(bidx)).values(sukeys))];

                [normmat,C,S]=normalize(rawmat,1,"range");
                %              [coeff,~,~]=pca(normmat);

                cv=cvpartition(min_trl,'KFold',10);
                dmat=struct();
                for trlType=lbls
                    dmat.(trlType).("band"+num2str(bidx))=cell2mat(cellfun(@(x) datasample(frmap.(trlType).("band"+num2str(bidx))(x),min_trl),sukeys,'UniformOutput',false));
                end

                for kf=1:cv.NumTestSets
                    [xx_train,yy_train,xx_test,yy_test]=deal([]);
                    for trlType=lbls
                        xx_train=[xx_train;dmat.(trlType).("band"+num2str(bidx))(training(cv,kf),:)];
                        xx_test=[xx_test;dmat.(trlType).("band"+num2str(bidx))(test(cv,kf),:)];
                        yy_train=[yy_train;repmat(trlType,nnz(training(cv,kf)),1)];
                        yy_test=[yy_test;repmat(trlType,nnz(test(cv,kf)),1)];
                    end
                    norm_train=normalize(xx_train,'center',C,'scale',S);
                    SVMM=fitcecoc(norm_train,yy_train);
                    norm_test=normalize(xx_test,'center',C,'scale',S);
                    modelPredict=SVMM.predict(norm_test);
                    result=[result;modelPredict==yy_test];
                end
            end
            out.(dec_tag).(sprintf('result_%dsu',n_su)).("band"+num2str(bidx))=result;
        end
    end
    save(sprintf('5_band_decoding_%s.mat',dec_tag),'out');
elseif opt.plot_dec
    load(sprintf('5_band_decoding_%s.mat',dec_tag),'out');
end

if opt.plot_dec
    switch opt.cmap
        case 'parula'
            cmap=flip(colormap(parula(5)));
        case 'cool'
            cmap=colormap(cool(5));
        otherwise
            cmap=colormap(turbo(5));
    end
    
    fh=figure('Color','w','Position',[100,100,720,450]);
    hold on;
    for bidx=1:(numel(ebound)-1)
        n_su=opt.n_su;
        phat=nan(1,numel(n_su));
        pci=nan(2,numel(n_su));
        for nidx=1:numel(n_su)
            dd=out.(dec_tag).(sprintf('result_%dsu',n_su(nidx))).("band"+num2str(bidx));
            [phat(nidx),pci(:,nidx)]=binofit(nnz(dd),numel(dd));
        end
        fill([n_su,fliplr(n_su)],[pci(1,:),fliplr(pci(2,:))],cmap(bidx,:),'EdgeColor','none','FaceAlpha',0.1);
        ph(bidx)=plot(n_su,phat,'Color',cmap(bidx,:),'LineStyle','-');
    end
    yline(0.5,'k--');
    xlim([0,500])
    su_tags=arrayfun(@(x) "Top "+num2str(x-20)+"%~"+num2str(x)+"%" ,20:20:100);
    legend(ph,su_tags,'Location','eastoutside','Orientation','vertical')
    xlabel('Number of neurons')
    ylabel('Classification accuracy (%)')
    ylim([0.25,1])
    set(gca(),'YTick',0:0.25:1,'YTickLabel',0:25:100)
    %     xlim([min(n_su),max(n_su)])
    title(dec_tag);
    %     exportgraphics(fh,sprintf('pct_decoding_%s.pdf',dec_tag),'ContentType','vector');
    exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);

end


end

% minp=min([sens_meta.wrs_p_d3(:,1:3),sens_meta.wrs_p_d6(:,1:3)],[],2);
% pts=prctile(minp,[10:10:90]);
% wave.pct_decoding(minp,pts,'n_su',[10,50,100,200,300,500],'lblidx',5)





%%
function relative_frac_heatmap()
sens_p=min([sens_meta.wrs_p_d3(:,1:3),sens_meta.wrs_p_d6(:,1:3)],[],2);
sens_p_win=[0,prctile(sens_p,[10:10:90]),1];

dur_p=min([dur_meta.wrs_p_s1(:,1:3),dur_meta.wrs_p_s2(:,1:3)],[],2);
dur_p_win=[0,prctile(dur_p,[10:10:90]),1];

count_mat=nan(10,10);
for s_idx=1:10
    for d_idx=1:10
        count_mat(s_idx,d_idx)=...
            nnz(sens_p>sens_p_win(s_idx) & ...
            sens_p<=sens_p_win(s_idx+1) & ...
            dur_p>dur_p_win(d_idx) & ...
            dur_p<=dur_p_win(d_idx+1));
    end
end

frac_mat=count_mat./sum(count_mat,'all');

fh=figure('Color','w');
imagesc(frac_mat,[0.002,0.02])
colormap('turbo')
cbh=colorbar();
set(gca(),'YDir','normal','XTick',0.5:1:10.5,'XTickLabel',0:10:100,...
    'YTick',0.5:1:10.5,'YTickLabel',0:10:100);
cbh.Label.String='Proportion of total population (%)';
cbh.TickLabels=cbh.Ticks*100;
xlabel('Olfactory coding rank, lower is better(%)');
ylabel('Duration coding rank, lower is better (%)');
grid on
grid minor
set(fh,'Position',[100,100,480,285])

end

function follow_up_stats()





eff_meta=ephys.effect_size_meta();


%% p-distribution histogram

sens_efsz=max(abs(eff_meta.cohen_d_olf),[],2);
sens_win=[min(sens_efsz)./2,prctile(sens_efsz,[20:20:100])];

dur_efsz=max(abs(eff_meta.cohen_d_dur),[],2);
dur_win=[min(dur_efsz)./2,prctile(dur_efsz,[20:20:100])];

cmap=flip(colormap(parula(5)));
bh=ephys.cmap_histogram(sens_efsz,sens_win,cmap);
set(gca(),'YScale','linear','XScale','linear')
set(gcf(),'Color','w')
set(gcf(),'Position',[100,100,520,100])
xlabel('p-value (permutation-test)')
ylabel('Number of neurons')
title('Olfactory guided 5 band binning')
% xlim([min(sens_efsz),max(sens_efsz)])
xlim([0,0.6])
exportgraphics(bh,'pct_decoding.pdf','ContentType','vector','Append',true);

cmap=colormap(cool(5));
bh=ephys.cmap_histogram(dur_efsz,dur_win,cmap);
set(gca(),'YScale','linear','XScale','linear')
set(gcf(),'Color','w')
set(gcf(),'Position',[100,100,520,100])
xlabel('p-value (permutation-test)')
ylabel('Number of neurons')
title('Duration guided 5 band binning')
xlim([0,0.6])
exportgraphics(bh,'pct_decoding.pdf','ContentType','vector','Append',true);

%% svm decoding
[fh,olf_dec_olf]=wave.pct_decoding(sens_efsz,sens_win,'n_su',[10,50,100,200,300,500],'lblidx',5,'cmap','parula','new_data',true,'calc_dec',true,'rpt',100)
[fh,dur_dec_dur]=wave.pct_decoding(dur_efsz,dur_win,'n_su',[10,50,100,200,300,500],'lblidx',8,'cmap','cool','new_data',true,'calc_dec',true,'rpt',100);

%% cross decoding
wave.pct_decoding(sens_efsz,sens_win,'n_su',[10,50,100,200,300,500],'lblidx',8,'cmap','parula','cross',true,'new_data',true,'calc_dec',true)
title('Rank by odor-encoding, decoding duration')
exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);
wave.pct_decoding(dur_efsz,dur_win,'n_su',[10,50,100,200,300,500],'lblidx',5,'cmap','cool','cross',true,'new_data',true,'calc_dec',true)
title('Rank by duration-encoding, decoding odor')
exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);



%% mixed coding heatmap
win_cnt=numel(sens_win)-1;
count_mat=nan(win_cnt,win_cnt);
for s_idx=1:win_cnt
    for d_idx=1:win_cnt
        count_mat(s_idx,d_idx)=...
            nnz(sens_efsz>sens_win(s_idx) & ...
            sens_efsz<=sens_win(s_idx+1) & ...
            dur_efsz>dur_win(d_idx) & ...
            dur_efsz<=dur_win(d_idx+1));
    end
end

frac_mat=count_mat./sum(count_mat,'all');

figure('Color','w')
imagesc(frac_mat,[0.02,0.06])
colormap('turbo')
cbh=colorbar();
set(gca(),'YDir','normal','XTick',0.5:1:10.5,'XTickLabel',0:20:100,...
    'YTick',0.5:1:10.5,'YTickLabel',0:20:100);
cbh.Label.String='Proportion of total population (%)';
cbh.TickLabels=cbh.Ticks*100;
xlabel('Olfactory coding rank, higher is better(%)');
ylabel('Duration coding rank, higher is better (%)');
exportgraphics(gcf(),'pct_decoding.pdf','ContentType','vector','Append',true);


%% Functional coupling
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
for bw=1:2
    pct_meta=struct();
    pct_meta.wave_id=zeros(size(sens_efsz));
    pct_meta.wave_id(sens_efsz<sens_win(1+bw) & dur_efsz<dur_win(1+bw))=1;
    pct_meta.wave_id(sens_efsz>=sens_win(6-bw) & dur_efsz<dur_win(1+bw))=2;
    pct_meta.wave_id(sens_efsz<sens_win(1+bw) & dur_efsz>=dur_win(6-bw))=3;
    pct_meta.wave_id(sens_efsz>=sens_win(6-bw) & dur_efsz>=dur_win(6-bw))=4;
    if false
        [fh,~,th]=bz.inter_wave_pct(pct_meta);
        title(th,sprintf('class of %d bins',bw));
        exportgraphics(fh,'pct_decoding.pdf','ContentType','vector','Append',true);
    end
    %% region-dist
    if true
        pct_map_cells=cell(1,4);
        grey_regs=ephys.getGreyRegs('range','grey');
        for mii=1:4
            pct_map_cells{mii}=containers.Map('KeyType','char','ValueType','any');
            for rr=grey_regs
                reg_sel=strcmp(meta.reg_tree(5,:),rr).';
                pct_map_cells{mii}(rr{1})=nnz(reg_sel & pct_meta.wave_id==mii)./nnz(reg_sel);
            end
        end
        fh=ephys.pct_proportion_GLM(pct_map_cells,'PearsonLogLog',...
            'range','grey','data_type',sprintf('class-%d bins',bw),'stats_type','pct',...
            'feat_tag',{'Class#1','Class#2','Class#3','Class#4'});
        for ff=reshape(fieldnames(fh),1,[])
            exportgraphics(fh.(ff{1}),'pct_decoding.pdf','ContentType','vector','Append',true);
        end
    end
end



end


