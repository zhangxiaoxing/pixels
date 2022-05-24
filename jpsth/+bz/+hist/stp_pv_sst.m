function fh=stp_pv_sst(stats,opt)
arguments
    stats (1,1) struct
    opt.binw (1,1) double = 200
    opt.dir (1,:) char {mustBeMember(opt.dir,{'from','to'})}  = 'from'
    opt.treedepth (1,:) double {mustBeInteger,mustBeInRange(opt.treedepth,1,7)} = 3:5;
    opt.relation (1,:) char {mustBeMember(opt.relation,{'same','diff'})} = 'diff'
    opt.title (1,:) char = []
    opt.minpair (1,1) double = 30
end

[~,~,ratiomap]=ref.get_pv_sst();
fh=figure('Color','w','Position',[100,100,400,300]);
hold on;
dcolors={'m','b','k'};
dsize=[16,49,100];
for currdepth=opt.treedepth
    if strcmp(opt.relation,'diff')
        congsel=stats.congru.maxiter(:,1)==0 & stats.congru.diff_reg(:,currdepth);
        nonmsel=stats.nonmem.maxiter(:,1)==0 & stats.nonmem.diff_reg(:,currdepth);
        incosel=stats.incongru.maxiter(:,1)==0 & stats.incongru.diff_reg(:,currdepth);
    else
        congsel=stats.congru.maxiter(:,1)==0 & stats.congru.same_reg(:,currdepth);
        nonmsel=stats.nonmem.maxiter(:,1)==0 & stats.nonmem.same_reg(:,currdepth);
        incosel=stats.incongru.maxiter(:,1)==0 & stats.incongru.same_reg(:,currdepth);
    end
    if strcmp(opt.dir,'from'), pos=1;else, pos=2;end
    congreg=cellfun(@(x) x{pos}{currdepth},stats.congru.reg(congsel,:),'UniformOutput',false);
    nonmreg=cellfun(@(x) x{pos}{currdepth},stats.nonmem.reg(nonmsel,:),'UniformOutput',false);
    incoreg=cellfun(@(x) x{pos}{currdepth},stats.incongru.reg(incosel,:),'UniformOutput',false);
    [GCcon,GRcon]=groupcounts(congreg);
    [GCnon,GRnon]=groupcounts(nonmreg);
    [GCinc,GRinc]=groupcounts(incoreg);
    
    % figure();histogram(GC,1:10:max(GC));
    congrudata=fliplr(stats.congru.postspk(congsel,2:end))*100;
    nonmdata=fliplr(stats.nonmem.postspk(nonmsel,2:end))*100;
    incodata=fliplr(stats.incongru.postspk(incosel,2:end))*100;
    
    con_gr=GRcon(GCcon>=opt.minpair);
    non_gr=GRnon(GCnon>=opt.minpair);
    inc_gr=GRinc(GCinc>=opt.minpair);
    
    % sel_gr=intersect(intersect(con_gr,non_gr),inc_gr);
    sel_gr=intersect(con_gr,non_gr);
    
    dim=ceil(sqrt(numel(sel_gr)));
    pairdata=[];
    for i=1:size(sel_gr)
        [~,congstats]=bootci(1000,@(x) sum(mean(x)),congrudata(strcmp(congreg,(sel_gr(i))),:));
        [~,nonmstats]=bootci(1000,@(x) sum(mean(x)),nonmdata(strcmp(nonmreg,(sel_gr(i))),:));
        [~,incostats]=bootci(1000,@(x) sum(mean(x)),incodata(strcmp(incoreg,(sel_gr(i))),:));
        diffc_n=congstats-nonmstats;
        mm=mean(diffc_n);
%         ci=bootci(1000,@(x) mean(x),diffc_n);
        scatter(ratiomap(sel_gr{i}),mm,dsize(6-currdepth),dcolors{currdepth-2},'filled','MarkerFaceAlpha',0.5);
        text(ratiomap(sel_gr{i}),mm,sel_gr{i},'VerticalAlignment','top','HorizontalAlignment','center');
        pairdata=vertcat(pairdata,[ratiomap(sel_gr{i}),mm]);
%         errorbar(i,mm,ci(1)-mm,ci(2)-mm,'k.');
        %         ylim([-1,10]);
        %         text(max(xlim()),max(ylim()),...
        %             sprintf('n= %d, %d, %d',nnz(strcmp(congreg,(sel_gr(i)))),nnz(strcmp(incoreg,(sel_gr(i)))),nnz(strcmp(nonmreg,(sel_gr(i))))),...
        %             'VerticalAlignment','top','HorizontalAlignment','right')
    end
%     set(gca(),'XTick',1:size(sel_gr),'XTickLabel',sel_gr,'XTickLabelRotation',90);
    if ~isempty(opt.title)
        sgtitle(opt.title)
    end
    [R,P]=corrcoef(pairdata(:,1),pairdata(:,2));
    text(max(xlim()),max(ylim()),sprintf('r = %.3f\np = %.3f',R(2),P(2)),'HorizontalAlignment','right','VerticalAlignment','top');
    
end