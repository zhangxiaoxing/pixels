% TODO migrate to wrs_mux data set
% major revision 20221104
function stats=sust_trans_correct_error(sel_meta,opt)
arguments
    sel_meta
    opt.single_bin (1,1) logical = true
    opt.plot_scatter (1,1) logical = false
    opt.plot_per_su (1,1) logical = true
    opt.plot_per_trial (1,1) logical = false
    opt.plot_showcase (1,1) logical = false
end
set(groot,'defaultTextFontSize',10);
sumeta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
stats=struct();
if opt.plot_scatter
    figure('Color','w','Position',[100,100,400,400]);
    hold on;
end

% should be selective types, i.e. olf, dur, mux
types={1:4,5:6,7:8};
type_desc={'mix','olf','dur'};
pref_trls={{'s1d3'};{'s1d6'};{'s2d3'};{'s2d3'};{'s1d3','s1d6'};{'s2d3','s2d6'};{'s1d3','s2d3'};{'s1d6','s2d6'}};

for typeIdx=1:numel(types)

    coding_sel=ismember(sel_meta.wave_id,types{typeIdx});
    usess=unique(sumeta.sess(coding_sel));

    % sust | transient selection previous version

    homedir=ephys.util.getHomedir('type','raw');
    [stats.(type_desc{typeIdx}).sust,stats.(type_desc{typeIdx}).transient]=deal(nan(0,4));
%     stats.(sprintf('type%d_pertrial',onetype))=cell(0,4); %skip for now
    switch typeIdx
        case 1
            sust_sel=(all(sel_meta.p_mux<0.05,2) | all(sel_meta.p_olf<0.05 | sel_meta.p_dur<0.05,2))...
                & ismember(sel_meta.wave_id,1:4);
            sig_bins=sel_meta.p_mux<0.05 | sel_meta.p_olf<0.05 | sel_meta.p_dur<0.05;
        case 2
            sust_sel=all(sel_meta.p_olf<0.05,2) & ismember(sel_meta.wave_id,5:6);
            sig_bins=sel_meta.p_olf<0.05;
        case 3
            sust_sel=all(sel_meta.p_olf<0.05,2) & ismember(sel_meta.wave_id,7:8);
            sig_bins=sel_meta.p_dur<0.05;
    end
    
    for sessIdx=reshape(usess,1,[])
        dpath=ephys.sessid2path(sessIdx);
        fpath=fullfile(homedir,dpath,'FR_All_1000.hdf5');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        fr=h5read(fpath,'/FR_All');
        
        trl.es1d3=find(trials(:,5)==4 & trials(:,8)==3 & trials(:,10)==0);
        trl.es2d3=find(trials(:,5)==8 & trials(:,8)==3 & trials(:,10)==0);
        trl.es1d6=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,10)==0);
        trl.es2d6=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,10)==0);

        trl.cs1d3=find(trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
        trl.cs2d3=find(trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
        trl.cs1d6=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);
        trl.cs2d6=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);

        if min([numel(trl.cs1d3);numel(trl.cs2d3);numel(trl.cs1d6);numel(trl.cs2d6);...
                numel(trl.es1d3);numel(trl.es2d3);numel(trl.es1d6);numel(trl.es2d6)],[],'all')<2
                continue
        end

        sess_sel=sumeta.sess==sessIdx;
        sessuid=sumeta.allcid(sess_sel & coding_sel);
        waveids=sel_meta.wave_id(sess_sel);
        sess_sust=sust_sel(sess_sel);
        sess_sig_bins=sig_bins(sess_sel,:);
        [~,su_pos]=ismember(sessuid,suid);
        for suidx=reshape(su_pos,1,[])
            su_waveid=waveids(suidx);
            su_pref_trls=pref_trls{su_waveid};
            if sess_sust(su_pos)
                bins=5:7;
            else
                bins=find(sess_sig_bins(suidx,:))+4;
            end

            [cpref,cnonpref,epref,enonpref]=deal([]);
            for trlfn=reshape(fieldnames(trl),1,[])
                % mean fr in selective bins

                trlType=char(trlfn);
                if startsWith(trlType,'c')
                    if ismember(replace(trlType,'c',''),su_pref_trls) %preferred correct
                        cpref=[cpref;reshape(fr(trl.(trlType),suidx,bins),[],1)];
                    else
                        cnonpref=[cnonpref;reshape(fr(trl.(trlType),suidx,bins),[],1)];
                    end
                else
                    if ismember(replace(trlType,'e',''),su_pref_trls) %preferred error
                        epref=[epref;reshape(fr(trl.(trlType),suidx,bins),[],1)];
                    else %nonpreferred error
                        enonpref=[enonpref;reshape(fr(trl.(trlType),suidx,bins),[],1)];
                    end
                end
            end

            delaymm=mean([mean(cpref),mean(cnonpref)]);
            delaystd=std([cpref;cnonpref]);
            if delaystd==0, continue;  end

            zcpref=(mean(cpref)-delaymm)./delaystd;
            zcnonpref=(mean(cnonpref)-delaymm)./delaystd; % for AUC curve maybe not necessary
            zepref=(mean(epref)-delaymm)./delaystd;
            zenonpref=(mean(enonpref)-delaymm)./delaystd;
            if sess_sust(suidx)
                stats.(type_desc{typeIdx}).sust=[stats.(type_desc{typeIdx}).sust;zcpref,zepref,zcnonpref,zenonpref];
            else
                stats.(type_desc{typeIdx}).transient=[stats.(type_desc{typeIdx}).transient;zcpref,zepref,zcnonpref,zenonpref];
            end

%             % skip for now
%             stats.(sprintf('type%d_pertrial',onetype))(idx,:)={cs1,cs2,es1,es2};


            %         if onetype==2 && mean(cs1)<mean(cs2),keyboard;end
            if opt.plot_showcase % perfcurve
                if min(cellfun(@(x) numel(x),{cs1,cs2,es1,es2}))<10
                    continue
                end
                [~,~,~,AUCC]=perfcurve([zeros(size(cs1));ones(size(cs2))],[cs1;cs2],0);
                [~,~,~,AUCE]=perfcurve([zeros(size(es1));ones(size(es2))],[es1;es2],0); 

                if abs(AUCC-0.5)>0.30 && abs(AUCE-0.5)<0.15
                    figure()
                    subplot(1,3,1)
                    hold on
                    cs1c=squeeze(mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6,suid==sumeta.allcid(ii),3:10),[1,2]));
                    cs2c=squeeze(mean(fr(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6,suid==sumeta.allcid(ii),3:10),[1 2]));
                    es1c=squeeze(mean(fr(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==6,suid==sumeta.allcid(ii),3:10),[1 2]));
                    es2c=squeeze(mean(fr(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==6,suid==sumeta.allcid(ii),3:10),[1 2]));

                    plot(cs1c,'-r')
                    plot(cs2c,'-b')
                    plot(es1c,':r')
                    plot(es2c,':b')
                    xline(1.5,'--k')


                    subplot(1,3,2)
                    hold on
                    histogram(cs1,-3:0.3:3,'FaceColor','r','FaceAlpha',0.4)
                    histogram(cs2,-3:0.3:3,'FaceColor','b','FaceAlpha',0.4)
                    subplot(1,3,3)
                    hold on
                    histogram(es1,-3:0.3:3,'FaceColor','r','FaceAlpha',0.4)
                    histogram(es2,-3:0.3:3,'FaceColor','b','FaceAlpha',0.4)
                    sgtitle(ii);
                    keyboard()
                end
            end

            %         if opt.plot_scatter
            %             ch(fidx)=scatter(mean(cs1),mean(cs2),36,colors{fidx},'MarkerFaceColor',colors{fidx},'MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');
            %             eh=scatter(mean(es1),mean(es2),16,'k','MarkerFaceColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
            %             plot([mean(cs1),mean(es1)],[mean(cs2),mean(es2)],':','LineWidth',0.5,'Color',colors{fidx});
            %         end


        end % per su
    end % per sessend
end

if opt.plot_scatter
    plot([-10,10],[-10,10],'--k')
    xlim([-2,2]);
    ylim(xlim());
    xlabel('Sample 1 normalized FR (Z-score)')
    ylabel('Sample 2 normalized FR (Z-score)')
    legend([ch(1),ch(2),eh],{'S1 correct trials','S2 correct trials','Error trials'});
end

% SU per session mean
% correct trial
if opt.plot_per_su
    for ff=["olf","mix","dur"]
        keyboard()
        stats.(ff).sust=stats.(ff).sust(stats.(ff).sust(:,1)>0,:);
        stats.(ff).transient=stats.(ff).transient(stats.(ff).transient(:,1)>0,:);

        fh=figure('Color','w','Position',[100,100,750,225]);
        tiledlayout(1,3)
        nexttile(1)
        hold on;
        ch=histogram(stats.(ff).sust(:,1),-0.7:0.1:1,'Normalization','probability','FaceColor','r','FaceAlpha',0.4);
        eh=histogram(stats.(ff).sust(:,2),-0.7:0.1:1,'Normalization','probability','FaceColor','k','FaceAlpha',0.4);
        xline(mean(stats.(ff).sust(:,1)),'--r','LineWidth',1);
        xline(mean(stats.(ff).sust(:,2)),'--k','LineWidth',1);
        xlim([-0.7,1]);
        title('Sustained')
        xlabel('Normalized FR (Z-Score)');
        ylabel('Probability')
        text(max(xlim()),max(ylim()),num2str(size(stats.(ff).sust,1)),'HorizontalAlignment','right','VerticalAlignment','top');
        [~,p]=ttest2(stats.(ff).sust(:,1),stats.(ff).sust(:,2));
        text(min(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','left','VerticalAlignment','top');

        nexttile(2)
        hold on;
        ch=histogram(stats.(ff).transient(:,1),-0.7:0.1:1,'Normalization','probability','FaceColor','r','FaceAlpha',0.4);
        eh=histogram(stats.(ff).transient(:,2),-0.7:0.1:1,'Normalization','probability','FaceColor','k','FaceAlpha',0.4);
        xline(mean(stats.(ff).transient(:,1)),'--r','LineWidth',1);
        xline(mean(stats.(ff).transient(:,2)),'--k','LineWidth',1);
        xlim([-0.7,1]);
        title('Transient')
        xlabel('Normalized FR (Z-Score)');
        ylabel('Probability')
        text(max(xlim()),max(ylim()),num2str(size(stats.(ff).transient,1)),'HorizontalAlignment','right','VerticalAlignment','top');
        [~,p]=ttest2(stats.(ff).transient(:,1),stats.(ff).transient(:,2));
        text(min(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','left','VerticalAlignment','top');

        %auc
%         onecolumn=size(stats.(sprintf('type%d',types(1))),1)+size(stats.(sprintf('type%d',types(2))),1);
%         [xc,yc,~,aucc]=perfcurve((1:2*onecolumn)>onecolumn,[stats.(sprintf('type%d',types(1)))(:,1);stats.(sprintf('type%d',types(2)))(:,2);stats.(sprintf('type%d',types(1)))(:,2);stats.(sprintf('type%d',types(2)))(:,1)],0);
%         [xe,ye,~,auce]=perfcurve((1:2*onecolumn)>onecolumn,[stats.(sprintf('type%d',types(1)))(:,3);stats.(sprintf('type%d',types(2)))(:,4);stats.(sprintf('type%d',types(1)))(:,4);stats.(sprintf('type%d',types(2)))(:,3)],0);
%         subplot(1,3,3);
%         hold on;
%         hc=plot(xc,yc,'-r','LineWidth',1);
%         he=plot(xe,ye,'-k','LineWidth',1);
%         legend([hc,he],{sprintf('Correct trials AUC=%0.3f',aucc),...
%             sprintf('Error trials AUC=%0.3f',auce)},...
%             'Location','southeast');
%         xlabel('False positive rate (fpr)');
%         ylabel('True positive rate (tpr)');
    end
%     if opt.single_bin
%         sgtitle(sprintf('%s averaged cross-trial bin %d',type,opt.bin));
%     else
%         sgtitle(sprintf('%s averaged cross-trial',type));
%     end
%     exportgraphics(fh,sprintf('%s_cross_trial_bin%d.pdf',type,opt.bin),'ContentType','vector');
end

%Per trial
if opt.plot_per_trial
    %correct trials
    figure('Color','w','Position',[100,100,1200,300]);
    subplot(1,3,1);
    hold on;
    xbins=-3:0.15:3;
    prefmat=cell2mat(cellfun(@(x) histcounts(x,xbins,'Normalization','probability'),...
        [stats.(sprintf('type%d_pertrial',types(1)))(:,1);...
        stats.(sprintf('type%d_pertrial',types(2)))(:,2)],...
        'UniformOutput',false));
    nonprefmat=cell2mat(cellfun(@(x) histcounts(x,xbins,'Normalization','probability'),...
        [stats.(sprintf('type%d_pertrial',types(1)))(:,2);...
        stats.(sprintf('type%d_pertrial',types(2)))(:,1)],...
        'UniformOutput',false));

    ph=bar(xbins(1:end-1)+0.075,mean(prefmat),'FaceColor','r','FaceAlpha',0.4);
    nph=bar(xbins(1:end-1)+0.075,mean(nonprefmat),'FaceColor','k','FaceAlpha',0.4);

    pcom=sum((xbins(1:end-1)+0.075).*mean(prefmat))./sum(mean(prefmat));
    npcom=sum((xbins(1:end-1)+0.075).*mean(nonprefmat))./sum(mean(nonprefmat));
    xline(pcom,'--r','LineWidth',1);
    xline(npcom,'--k','LineWidth',1);
    title('Correct trials')
    xlabel('Normalized FR (Z-Score)');
    ylabel('Probability')

    %error trials
    subplot(1,3,2);
    hold on;
    xbins=-3:0.15:3;
    prefmat=cell2mat(cellfun(@(x) histcounts(x,xbins,'Normalization','probability'),...
        [stats.(sprintf('type%d_pertrial',types(1)))(:,3);...
        stats.(sprintf('type%d_pertrial',types(2)))(:,4)],...
        'UniformOutput',false));
    nonprefmat=cell2mat(cellfun(@(x) histcounts(x,xbins,'Normalization','probability'),...
        [stats.(sprintf('type%d_pertrial',types(1)))(:,4);...
        stats.(sprintf('type%d_pertrial',types(2)))(:,3)],...
        'UniformOutput',false));

    ph=bar(xbins(1:end-1)+0.075,nanmean(prefmat),'FaceColor','r','FaceAlpha',0.4);
    nph=bar(xbins(1:end-1)+0.075,nanmean(nonprefmat),'FaceColor','k','FaceAlpha',0.4);

    pcom=sum((xbins(1:end-1)+0.075).*nanmean(prefmat))./sum(nanmean(prefmat));
    npcom=sum((xbins(1:end-1)+0.075).*nanmean(nonprefmat))./sum(nanmean(nonprefmat));
    xline(pcom,'--r','LineWidth',1);
    xline(npcom,'--k','LineWidth',1);

    title('Error trials');
    legend([ph,nph],{'Prefered','Non-prefered'});
    xlabel('Normalized FR (Z-Score)');
    ylabel('Probability')

    %auc

    subplot(1,3,3);
    hold on;
    xq=0.005:0.01:1;
    prefs1=arrayfun(@(x) auc_per_su(stats.(sprintf('type%d_pertrial',types(1)))(x,:),0,xq),...
        1:size(stats.(sprintf('type%d_pertrial',types(1))),1));
    prefs2=arrayfun(@(x) auc_per_su(stats.(sprintf('type%d_pertrial',types(2)))(x,:),1,xq),...
        1:size(stats.(sprintf('type%d_pertrial',types(2))),1));

    hc=plot(xq,mean([cell2mat({prefs1.yc}.');cell2mat({prefs2.yc}.')]),'-r','LineWidth',1);
    he=plot(xq,mean([cell2mat({prefs1.ye}.');cell2mat({prefs2.ye}.')]),'-k','LineWidth',1);
    legend([hc,he],{sprintf(...
        'Correct trials AUC=%0.2f',mean([prefs1.aucc,prefs2.aucc]))...
        ,sprintf(...
        'Error trials AUC=%0.2f',mean([prefs1.auce,prefs2.auce]))},...
        'Location','southeast');
    xlabel('False positive rate (fpr)');
    ylabel('True positive rate (tpr)');
    sgtitle(sprintf('%s stats from trials, averaged over neurons',type));
end

end


function out=auc_per_su(data,pos,xq)
out=struct();
labels=(1:numel([data{1};data{2}]))>numel(data{1});
scores=[data{1};data{2}];
[xc,yc,~,out.aucc]=perfcurve(labels,scores,pos);
[G,uxc]=findgroups(xc);
mmy=splitapply(@(x) mean(x), yc, G);
out.yc=interp1(uxc,mmy,xq);

labels=(1:numel([data{3};data{4}]))>numel(data{3});
scores=[data{3};data{4}];
[xe,ye,~,out.auce]=perfcurve(labels,scores,pos);
[G,uxe]=findgroups(xe);
mmy=splitapply(@(x) mean(x), ye, G);
out.ye=interp1(uxe,mmy,xq);

end
