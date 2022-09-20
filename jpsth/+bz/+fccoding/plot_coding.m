function plot_coding(sel_meta,opt)
arguments
    sel_meta
    opt.plot_trial_frac (1,1) logical = false
    opt.plot_fwd_rev (1,1) logical = false
    opt.plot_coding_idx (1,1) logical = false
    opt.plot_coding_idx_shuf (1,1) logical = false
    opt.plot_svm (1,1) logical = true
    opt.wave_dir (1,1) logical = false
end

persistent metas stats per_trial_norm

% [~,~,ratiomap]=ref.get_pv_sst();
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
load('OBM1Map.mat','OBM1map')

% {FWD/RWD,H2L/L2H/Local, S1/S2*correct/error*3s/6s}
% congru, incongru, nonmem

if opt.plot_trial_frac
    [metas,stats,~]=bz.fccoding.get_fc_coding(pct_meta,'no_jitter',false); 
    % meta:sess,suid1,suid2,reg:2;
    %TODO assign mem_type_selection 
    error('Needs update');
    congrus1=ismember(metas(:,4),1:2) & ismember(metas(:,5),1:2) & all(~ismissing(stats),2);
    congrus2=ismember(metas(:,4),3:4) & ismember(metas(:,5),3:4) & all(~ismissing(stats),2);
    ci=bootci(500,@(x) mean(x),[stats(congrus1,5);stats(congrus2,6)]);
    mm=mean([stats(congrus1,5);stats(congrus2,6)]);
    fh=figure('Color','w','Position',[100,100,100,235]);
    hold on
    swarmchart(ones(nnz(congrus1+congrus2),1),[stats(congrus1,5);stats(congrus2,6)],1,[0.8,0.8,0.8],'o','MarkerFaceAlpha',0.2,'MarkerFaceColor','k','MarkerEdgeColor','none')
    boxplot([stats(congrus1,5);stats(congrus2,6)],'Colors','k','Symbol','','Widths',0.8)
%     errorbar(1,mm,ci(1)-mm,ci(2)-mm,'k.','CapSize',30);
    set(gca(),'XTick',[],'YTick',0:0.2:1,'YTickLabel',0:20:100)
    xlabel('Real-shuffle')
    ylabel('Appeared trials (%)')
    exportgraphics(fh,'FC_appearance.pdf')
end

if opt.plot_fwd_rev
    [metas,stats,fwd_rev]=bz.fccoding.get_fc_coding('no_jitter',false);
    congrus1=ismember(metas(:,4),1:2) & ismember(metas(:,5),1:2) & all(~ismissing(stats),2);
    congrus2=ismember(metas(:,4),3:4) & ismember(metas(:,5),3:4) & all(~ismissing(stats),2);
    fridx1=[fwd_rev(congrus1,1)];
    fridx2=[fwd_rev(congrus2,2)];
    fh=figure('Color','w','Position',[100,100,100,235]);
    hold on
    swarmchart(ones(nnz(congrus1+congrus2),1),[stats(congrus1,5);stats(congrus2,6)],4,[0.8,0.8,0.8],'o','MarkerFaceAlpha',0.2,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','none')
    boxplot([stats(congrus1,5);stats(congrus2,6)],'Colors','k','Symbol','','Widths',0.8)
    set(gca(),'XTick',[],'YTick',0:0.2:1,'YTickLabel',0:20:100)
    xlabel('Real-shuffle')
    ylabel('Appeared trials (%)')
    exportgraphics(fh,'FC_appearance.pdf')
end


if opt.plot_coding_idx
    [metas,stats,~]=bz.fccoding.get_fc_coding('no_jitter',true);
    out_same=plot_one('same',metas,stats,idmap,OBM1map,false);
    out_h2l=plot_one('h2l',metas,stats,idmap,OBM1map,false);
    out_l2h=plot_one('l2h',metas,stats,idmap,OBM1map,false);
    
    fh=figure('Color','w','Position',[100,100,235,235]);
    hold on;
    for pi=0:3:6
        if pi==0,mm=out_same.mm;cic=out_same.cic;cie=out_same.cie;...
        elseif pi==3,mm=out_l2h.mm;cic=out_l2h.cic;cie=out_l2h.cie;...
        else,mm=out_h2l.mm;cic=out_h2l.cic;cie=out_h2l.cie;
        end
        ch=bar(pi+1,mm(1),'FaceColor','w','EdgeColor','k');
        eh=bar(pi+2,mm(2),'FaceColor','k','EdgeColor','k');
        errorbar(pi+(1:2),mm,[cic(1),cie(1)]-mm,[cic(2),cie(2)]-mm,'k.');
    end
    ylabel('F.C. selectivity  index')
%     text(max(xlim()),max(ylim()),sprintf('p=%.3f',ranksum([dists1;dists2],[dists1e;dists2e])),'HorizontalAlignment','right','VerticalAlignment','top');
    ylim([0,0.25]);
    set(gca,'XTick',1.5:3:7.5,'XTickLabel',{'Within reg.','Sens. to motor.','Motor to sens.'},'XTickLabelRotation',45)
    legend([ch,eh],{'Correct trials','Error trials'});
    pc=anovan([out_same.datac;out_h2l.datac;out_l2h.datac],{[zeros(size(out_same.datac));ones(size(out_h2l.datac));2*ones(size(out_l2h.datac))]})
    pe=anovan([out_same.datae;out_h2l.datae;out_l2h.datae],{[zeros(size(out_same.datae));ones(size(out_h2l.datae));2*ones(size(out_l2h.datae))]})
    exportgraphics(fh,'FC_coding_hier.pdf')
    keyboard()
end

if opt.plot_svm
    if isempty(metas) ||isempty(per_trial_norm)
        [metas,~,~,per_trial]=bz.fccoding.get_fc_coding(sel_meta,'no_jitter',true,'per_trial',true);

        %normalize
        [centt,stdd]=deal(nan(size(per_trial(:,1))));
        per_trial_norm=cell(size(per_trial(:,1:2)));
        for ii=1:numel(centt)
            [N,centt(ii),stdd(ii)]=normalize([per_trial{ii,1};per_trial{ii,2}]);
            per_trial_norm{ii,1}=N(1:numel(per_trial{ii,1}));
            per_trial_norm{ii,2}=N(numel(per_trial{ii,1})+1:end);
            per_trial_norm{ii,3}=normalize(per_trial{ii,3},'center',centt(ii),'scale',stdd(ii));
            per_trial_norm{ii,4}=normalize(per_trial{ii,4},'center',centt(ii),'scale',stdd(ii));
        end
    end

    congrusel=all(ismember(metas(:,4:5),[1,5]),2) | all(ismember(metas(:,4:5),[3,5]),2) ...
        | all(ismember(metas(:,4:5),[2,6]),2) | all(ismember(metas(:,4:5),[4,6]),2);
    if opt.wave_dir

        repreg=cat(3,repmat(metas(:,6),1,6),repmat(metas(:,7),1,6));
        [~,sig_same,sig_h2l,sig_l2h]=bz.util.diff_at_level(repreg,'hierarchy',true,'hiermap','AON','descend',true);
        mtypes={congrusel & sig_l2h(:,5),...
            congrusel & sig_h2l(:,5),...
            congrusel & sig_same(:,5),...
            any(ismember(metas(:,4:5),[1,3,5]),2) & any(ismember(metas(:,4:5),[2,4,6]),2),...
            all(metas(:,4:5)==0,2)};

        nsu_grp=10:10:30;
    else
        mtypes={congrusel,...
            any(ismember(metas(:,4:5),[1,3,5]),2) & any(ismember(metas(:,4:5),[2,4,6]),2),...
            all(metas(:,4:5)==0,2)};
        nsu_grp=25:25:100;
    end


    NTRIAL=40;
    NERR=5;
    NRPT=50;
    dec_result=struct();
    
    for NSU=nsu_grp
        [decdata.s1,decdata.s2]=deal(nan(NSU,NTRIAL,1,NRPT));
        [decdata.e1,decdata.e2]=deal(nan(NSU,NERR,1,NRPT));

        result=cell(1,3);
        % congru, incong, nonmem
        % TODO wave direction

        for midx=1:numel(mtypes)
            su_sel=find(min(cellfun(@(x) size(x,1),per_trial_norm(:,1:2)),[],2)>NTRIAL & mtypes{midx} & centt>1 & stdd>0 & min(cellfun(@(x) size(x,1),per_trial(:,3:4)),[],2)>NERR);
            disp([midx,numel(su_sel)]);
            for rpt=1:NRPT
                su_perm=randsample(su_sel,NSU);
                decdata.s1(:,:,1,rpt)=cell2mat(arrayfun(@(x) randsample(per_trial_norm{x,1},NTRIAL),su_perm,'UniformOutput',false).').';
                decdata.s2(:,:,1,rpt)=cell2mat(arrayfun(@(x) randsample(per_trial_norm{x,2},NTRIAL),su_perm,'UniformOutput',false).').';
                decdata.e1(:,:,1,rpt)=cell2mat(arrayfun(@(x) randsample(per_trial_norm{x,3},NERR),su_perm,'UniformOutput',false).').';
                decdata.e2(:,:,1,rpt)=cell2mat(arrayfun(@(x) randsample(per_trial_norm{x,4},NERR),su_perm,'UniformOutput',false).').';
            end
            if midx==1
                result{midx}=fc.dec.dec(decdata,'decoder','SVM','svmcost',1);
            else
                result{midx}=fc.dec.dec(decdata,'decoder','SVM','svmcost',1);
            end
        end
        dec_result.(sprintf('SU%d',NSU))=result;
    end
    
    
    colors={'r','b','k','c','m'};
    fh=figure('Color','w','Position',[100,100,235,235]);
    hold on;
    for mi=1:1:numel(mtypes)
        cvmat=cell2mat(arrayfun(@(x) dec_result.(sprintf('SU%d',x)){mi}.cvcorr{1},nsu_grp,'UniformOutput',false));
        mm=mean(cvmat).*100;
        cim=bootci(500,@(x) mean(x).*100,cvmat);
        ph(mi)=plot(nsu_grp,mm,strjoin({'-',colors{mi}},''));
        fill([nsu_grp,fliplr(nsu_grp)],[cim(1,:),fliplr(cim(2,:))],colors{mi},'EdgeColor','none','FaceAlpha',0.1);
    end
    
    cvmat=cell2mat(arrayfun(@(x) dec_result.(sprintf('SU%d',x)){mi}.errcorr{1},nsu_grp,'UniformOutput',false));
    mm=mean(cvmat).*100;
    cim=bootci(500,@(x) mean(x).*100,cvmat);
    ph(numel(mtypes)+1)=plot(nsu_grp,mm,'--r');
    fill([nsu_grp,fliplr(nsu_grp)],[cim(1,:),fliplr(cim(2,:))],'r','EdgeColor','none','FaceAlpha',0.1);

    ylabel('Odor classification accuracy');
    xlabel('Number of FC pairs');
    %     xlim([50,100])
    if opt.wave_dir
        legend(ph,{'Upward','Downward','Same-region','Diff-mem','Nonmemory','(Error trial)'},Location="northoutside");
    else
        legend(ph,{'Same memory','Diff. memory','Nonmemory','(Error trial)'},Location="northoutside");
    end
%     exportgraphics(fh,'FC_coding_SVM.pdf')
end



if opt.plot_coding_idx_shuf
    [metas,stats,~]=bz.fccoding.get_fc_coding('no_jitter',true,'shuffle',false);
    out_any=plot_one('any',metas,stats,idmap,OBM1map,false);
    keyboard() 
    %% could be very time consuming for large repeats!
    %should backport from cluster.
    parfor rpt=1:500
        [metas,stats,~]=bz.fccoding.get_fc_coding('no_jitter',true,'shuffle',true);
        out_shuf(rpt)=plot_one('shuffle',metas,stats,idmap,OBM1map,false);
    end
    shufdatac=cell2mat({out_shuf.datac}.');
    shufmm=mean(shufdatac);
    shufci=bootci(500,@(x) mean(x), shufdatac);
    
    mm=[out_any.mm,shufmm];
    ci=[out_any.cic,out_any.cie,shufci];
    p=anovan([out_any.datac;out_any.datae;shufdatac],{[zeros(size(out_any.datac));ones(size(out_any.datae));2*ones(size(shufdatac))]})
    %% following import data
    fh=figure('Color','w','Position',[100,100,235,235]);
    hold on;
    bh=bar(mm,'FaceColor','w','EdgeColor','k');
    errorbar(bh.XEndPoints,mm,ci(1,:)-mm,ci(2,:)-mm,'k.');
    ylabel('F.C. spike selectivity  index')
%     text(max(xlim()),max(ylim()),sprintf('p=%.3f',ranksum([dists1;dists2],[dists1e;dists2e])),'HorizontalAlignment','right','VerticalAlignment','top');
    ylim([0,0.17]);
    set(gca,'XTick',1:3,'XTickLabel',{'Correct','Error','Shuffle'},'XTickLabelRotation',45)
%     legend([ch,eh],{'Correct trials','Error trials'});
%     pc=anovan([out_same.datac;out_h2l.datac;out_l2h.datac],{[zeros(size(out_same.datac));ones(size(out_h2l.datac));2*ones(size(out_l2h.datac))]})
%     pe=anovan([out_same.datae;out_h2l.datae;out_l2h.datae],{[zeros(size(out_same.datae));ones(size(out_h2l.datae));2*ones(size(out_l2h.datae))]})
    exportgraphics(fh,'FC_coding_hier.pdf')
    keyboard()
end


end
function out=plot_one(type,meta_in,stat_in,idmap,index_map,plot_)
w_key=arrayfun(@(x) index_map.isKey(char(idmap.ccfid2reg(x))),meta_in(:,6:7));
metas=meta_in(all(w_key,2),:); % only region in map kept below
stats=stat_in(all(w_key,2),:);
index_value=arrayfun(@(x) index_map(char(idmap.ccfid2reg(x))),metas(:,6:7));
switch type
    case 'same'
        congrus1=all(ismember(metas(:,4:5),1:2),2) & all(~ismissing(stats),2) & metas(:,6)==metas(:,7) ;
        congrus2=all(ismember(metas(:,4:5),3:4),2) & all(~ismissing(stats),2) & metas(:,6)==metas(:,7);
        ftitle='Within region';
    case 'l2h'
        congrus1=all(ismember(metas(:,4:5),1:2),2) & all(~ismissing(stats),2) & index_value(:,1)>index_value(:,2);
        congrus2=all(ismember(metas(:,4:5),3:4),2) & all(~ismissing(stats),2) & index_value(:,1)>index_value(:,2);
        ftitle='Lower to higher';
    case 'h2l'
        congrus1=all(ismember(metas(:,4:5),1:2),2) & all(~ismissing(stats),2) & index_value(:,1)<index_value(:,2);
        congrus2=all(ismember(metas(:,4:5),3:4),2) & all(~ismissing(stats),2) & index_value(:,1)<index_value(:,2);
        ftitle='Higher to lower';
    case 'any' % only region in map kept, see above
        congrus1=all(ismember(metas(:,4:5),1:2),2) & all(~ismissing(stats),2);
        congrus2=all(ismember(metas(:,4:5),3:4),2) & all(~ismissing(stats),2);
    case 'shuffle' % only region in map kept, see above
        congrus1=all(ismember(metas(:,4:5),1:2),2) & all(~ismissing(stats),2);
        congrus2=all(ismember(metas(:,4:5),3:4),2) & all(~ismissing(stats),2);        
end

if ~strcmp(type,'shuffle')
dists1=arrayfun(@(x) -diff(stats(x,1:2),1,2)./sum(stats(x,1:2),2),find(congrus1 & sum(stats(:,1:2),2)>10));
dists2=arrayfun(@(x) diff(stats(x,1:2),1,2)./sum(stats(x,1:2),2),find(congrus2 & sum(stats(:,1:2),2)>10));
dists1e=arrayfun(@(x) -diff(stats(x,3:4),1,2)./sum(stats(x,3:4),2),find(congrus1 & sum(stats(:,3:4),2)>10));
dists2e=arrayfun(@(x) diff(stats(x,3:4),1,2)./sum(stats(x,3:4),2),find(congrus2& sum(stats(:,3:4),2)>10));

mm=[mean([dists1;dists2]),mean([dists1e;dists2e])];
cic=bootci(1000,@(x) mean(x),[dists1;dists2]);
cie=bootci(1000,@(x) mean(x),[dists1e;dists2e]);
out.mm=mm;
out.cic=cic;
out.cie=cie;
out.wrsp=ranksum([dists1;dists2],[dists1e;dists2e]);
out.datac=[dists1;dists2];
out.datae=[dists1e;dists2e];
else
    dists1=arrayfun(@(x) -diff(stats(x,1:2),1,2)./sum(stats(x,1:2),2),find(congrus1 & sum(stats(:,1:2),2)>10));
    dists2=arrayfun(@(x) diff(stats(x,1:2),1,2)./sum(stats(x,1:2),2),find(congrus2 & sum(stats(:,1:2),2)>10));
    mm=mean([dists1;dists2]);
    out.mm=mm;
    out.datac=[dists1;dists2];
end

if plot_
  
    fh=figure('Color','w','Position',[100,100,900,300]);
    subplot(1,3,1);
    hold on;
    ph1=scatter(stats(congrus1,1),stats(congrus1,2),9,'MarkerFaceColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
    ph2=scatter(stats(congrus2,1),stats(congrus2,2),9,'MarkerFaceColor','b','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
    plot([-10,50],[-10,50],'--k')
    xlim([-10,50])
    ylim([-10,50])
    title('Correct trials')
    legend([ph1,ph2],{'S1 congruent','S2 congruent'},'Location','northoutside');
    xlabel('S1 F.C. real-shuffle');
    ylabel('S2 F.C. real-shuffle');
    subplot(1,3,2);
    hold on;
    ph1=scatter(stats(congrus1,3),stats(congrus1,4),9,'MarkerFaceColor','m','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
    ph2=scatter(stats(congrus2,3),stats(congrus2,4),9,'MarkerFaceColor','c','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
    plot([-10,50],[-10,50],'--k')
    xlim([-10,50])
    ylim([-10,50])
    title('Error trials')
    legend([ph1,ph2],{'S1 congruent','S2 congruent'},'Location','northoutside');
    xlabel('S1 F.C. real-shuffle');
    ylabel('S2 F.C. real-shuffle');
    
    
    subplot(1,3,3);
    hold on;
    % ch=histogram([dists1;dists2],-6:0.5:6,'Normalization','probability');
    % eh=histogram([dists1e;dists2e],-6:0.5:6,'Normalization','probability');
    
    bar(1,mm(1),'FaceColor','w','EdgeColor','k')
    bar(2,mm(2),'FaceColor','k','EdgeColor','k')
    errorbar(1:2,mm,[cic(1),cie(1)]-mm,[cic(2),cie(2)]-mm,'k.');
    set(gca(),'XTick',1:2,'XTickLabel',{'Correct','Error'},'XTickLabelRotation',45)
    ylabel('F.C. selectivity  index')
    text(max(xlim()),max(ylim()),sprintf('p=%.3f',ranksum([dists1;dists2],[dists1e;dists2e])),'HorizontalAlignment','right','VerticalAlignment','top');
    ylim([-0.1,0.7]);
    sgtitle(ftitle);
end

end

function out=rnd_trial(incell,ntrial)
su_perm=randperm(numel(incell{1},1));
out=incell{su_perm(1:ntrial)};
end
