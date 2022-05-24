function pct_mem_neuron(rings,opt)
% 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched
arguments
    rings (:,3) cell
    opt.subsel (1,1) logical = false
end

sess_tagged=arrayfun(@(x) {x*100000+rings{x,1},...
    x*100000+rings{x,2},...
    x*100000+rings{x,3}},...
    1:size(rings,1),'UniformOutput',false)';

types=cell(numel(sess_tagged),3);
for ss=1:size(sess_tagged,1)
    for sidx=1:3
        rr=sess_tagged{ss}{sidx};
        if ~isempty(rr)
            types{ss,sidx}=arrayfun(@(x) bz.rings.ucid2memtype(x),rr);
        end
    end
end

fh=figure('Color','w','Position',[100,100,540,180]);

for sidx=1:3
    onesizemat=cell2mat(types(cellfun(@(x) ~isempty(x),types(:,sidx)),sidx));
    onesizemat=onesizemat(all(onesizemat>=0,2),:);
    frac=sum(onesizemat>0,2)./size(onesizemat,2);
    [GC,GR]=groupcounts(frac);
    [~,bootsam]=bootstrp(100,[],frac);
    bootstat=cell2mat(arrayfun(@(bi) ...
        arrayfun(@(x) nnz(frac(bootsam(:,bi))==x),GR),1:size(bootsam,2),...
        'UniformOutput',false))';
    ci=prctile(bootstat,[2.5,97.5])./numel(frac)*100;
    mm=GC'./numel(frac)*100;
    subplot(1,3,sidx);
    hold on;
    bar(GR*100,mm,'FaceColor','none','EdgeColor','k','LineWidth',1);
    errorbar(GR*100,mm,ci(1,:)-mm,ci(2,:)-mm,'k.','LineWidth',1)
    ylim([0,50])
    ylabel('Normalized proportion (%)')
    xlabel('Fraction of memory-cell (%)');
    set(gca,'FontSize',10);
end
exportgraphics(fh,fullfile('bzdata','pct_mem_neuron.pdf'));
end