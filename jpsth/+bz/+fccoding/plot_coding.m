function plot_coding(sel_meta,opt)
arguments
    sel_meta
    opt.plot_trial_frac (1,1) logical = false
    opt.plot_fwd_rev (1,1) logical = false
    opt.plot_coding_idx (1,1) logical = false
    opt.plot_coding_idx_shuf (1,1) logical = false
    opt.plot_svm (1,1) logical = true
    opt.wave_dir (1,1) logical = false
    opt.dtype (1,:) char {mustBeMember(opt.dtype,{'olf','dur'})} ='olf'
    opt.nrpt (1,1) double = 100
end

% [~,~,ratiomap]=ref.get_pv_sst();
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
% load('OBM1Map.mat','OBM1map')

% {FWD/RWD,H2L/L2H/Local, S1/S2*correct/error*3s/6s}
% congru, incongru, nonmem

if opt.plot_svm
    dtype=opt.dtype;
    nfc_grp=[10 50 100 250 500 750 1000];

    NTRIAL=20;
    NERR=2;
    NRPT=opt.nrpt;
    dec_result=struct();
    % resample trials
    stats_all=cell(1,NRPT/5);
    for ii=1:NRPT/5
        [~,stats_all{ii}]=bz.fccoding.get_fc_coding(sel_meta,'force_update',true,'type',dtype); % olfactory duration mix
        [~,stats_incong{ii}]=bz.fccoding.get_fc_coding(sel_meta,'force_update',true,'type',dtype,'incong',true); % olfactory duration mix
        
    end
    S=max(cell2mat(cellfun(@(x) max([x.lbl1;x.lbl2]).',stats_all,'UniformOutput',false)),[],2).';
    S(S==0)=1;

    Sinc=max(cell2mat(cellfun(@(x) max([x.lbl1;x.lbl2]).',stats_incong,'UniformOutput',false)),[],2).';
    Sinc(Sinc==0)=1;


    for N_pair=nfc_grp
        [decdata.s1,decdata.s2,incongdata.s1,incongdata.s2]=deal(nan(N_pair,NTRIAL,1,NRPT));
        [decdata.e1,decdata.e2,incongdata.e1,incongdata.e2]=deal(nan(N_pair,NERR,1,NRPT));
%             result=cell(1,3);
        % congru, incong, nonmem
        % TODO wave direction

        % TODO multiple types of pairs, i.e. congru, in congru ,same
        % region, different region etc.
        %         for midx=1:numel(mtypes)
        for rpt=1:NRPT
            ridx=ceil(rpt/5);
            stats=stats_all{ridx};
            incong=stats_incong{ridx};
            for ff=reshape(fieldnames(stats),1,[])
                stats.(char(ff))=stats.(char(ff))./repmat(S,size(stats.(char(ff)),1),1);
                incong.(char(ff))=incong.(char(ff))./repmat(Sinc,size(incong.(char(ff)),1),1);
            end
            fc_perm=randsample(size(stats.lbl1,2),N_pair);
            decdata.s1(:,:,1,rpt)=stats.lbl1(:,fc_perm).';
            decdata.s2(:,:,1,rpt)=stats.lbl2(:,fc_perm).';
            decdata.e1(:,:,1,rpt)=stats.e1(:,fc_perm).';
            decdata.e2(:,:,1,rpt)=stats.e2(:,fc_perm).';

            fc_perm=randsample(size(incong.lbl1,2),N_pair);
            incongdata.s1(:,:,1,rpt)=incong.lbl1(:,fc_perm).';
            incongdata.s2(:,:,1,rpt)=incong.lbl2(:,fc_perm).';
            incongdata.e1(:,:,1,rpt)=incong.e1(:,fc_perm).';
            incongdata.e2(:,:,1,rpt)=incong.e2(:,fc_perm).';
        end
        result=fc.dec.dec(decdata,'decoder','SVM','svmcost',1);
        inc_result=fc.dec.dec(incongdata,'decoder','SVM','svmcost',1);
        %         end % TODO multiple pair types
        dec_result.(sprintf('FC%d',N_pair))=result;
        incong_result.(sprintf('FC%d',N_pair))=inc_result;
    end


    fh=figure('Color','w','Position',[100,100,235,235]);
    hold on;

    cvmat=cell2mat(arrayfun(@(x) dec_result.(sprintf('FC%d',x)).shufcorr{1},nfc_grp,'UniformOutput',false));
    [mm,cim]=binofit(sum(cvmat),repmat(size(cvmat,1),1,size(cvmat,2)));
    phe=plot(nfc_grp,mm,'-k');
    fill([nfc_grp,fliplr(nfc_grp)],[cim(:,1);flip(cim(:,2))],'k','EdgeColor','none','FaceAlpha',0.1);

    cvmat=cell2mat(arrayfun(@(x) dec_result.(sprintf('FC%d',x)).errcorr{1},nfc_grp,'UniformOutput',false));
    [mm,cim]=binofit(sum(cvmat),repmat(size(cvmat,1),1,size(cvmat,2)));
    phe=plot(nfc_grp,mm,'-b');
    fill([nfc_grp,fliplr(nfc_grp)],[cim(:,1);flip(cim(:,2))],'b','EdgeColor','none','FaceAlpha',0.1);
    
    cvmat=cell2mat(arrayfun(@(x) dec_result.(sprintf('FC%d',x)).cvcorr{1},nfc_grp,'UniformOutput',false));
    [mm,cim]=binofit(sum(cvmat),repmat(size(cvmat,1),1,size(cvmat,2)));
    phc=plot(nfc_grp,mm,'-r');
    fill([nfc_grp,fliplr(nfc_grp)],[cim(:,1);flip(cim(:,2))],'r','EdgeColor','none','FaceAlpha',0.1);

    % incongru
    cvmat=cell2mat(arrayfun(@(x) incong_result.(sprintf('FC%d',x)).cvcorr{1},nfc_grp,'UniformOutput',false));
    [mm,cim]=binofit(sum(cvmat),repmat(size(cvmat,1),1,size(cvmat,2)));
    phi=plot(nfc_grp,mm,'-m');
    fill([nfc_grp,fliplr(nfc_grp)],[cim(:,1);flip(cim(:,2))],'m','EdgeColor','none','FaceAlpha',0.1);


    ylabel('Classification accuracy');
    xlabel('Number of FC pairs');
    title(dtype)
    ylim([0.4,1])


end

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



if opt.plot_fwd_rev % probably broken as of 22.11.10
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

if opt.plot_coding_idx % TODO might work with minor updates
    [metas,stats]=bz.fccoding.get_fc_coding('no_jitter',true);
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
