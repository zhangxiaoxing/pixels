function fh=ring_dist(sig,pair,rings,msize,opt)
arguments
    sig (1,1) struct
    pair (1,1) struct
    rings (:,3) cell
    msize (1,1) double {mustBeMember(msize,3:5)}
    opt.type (1,:) char {mustBeMember(opt.type,{'max','min'})} = 'max'
end

addpath('k:\code\align\')
[sig.reg_dist,pair.reg_dist]=get_conn_dist(sig,pair);

dist_span=cell(size(rings,1),1);
for si=1:size(rings,1)
    if ~isempty(rings{si,msize-2})
        dist_one=nan(size(rings{si,msize-2},1),msize);
        for ri=1:size(rings{si,msize-2},1)
            one_ring=rings{si,msize-2}(ri,:);
            for di=1:(size(one_ring,2)-1)
                dist_one(ri,di)=sig.reg_dist(sig.sess==si & sig.suid(:,1)==one_ring(di) & sig.suid(:,2)==one_ring(di+1));
            end
            dist_one(ri,end)=sig.reg_dist(sig.sess==si & sig.suid(:,1)==one_ring(end) & sig.suid(:,2)==one_ring(1));
        end
    else
        dist_one=[];
    end
    dist_span{si}=dist_one;
end
span_mat=cell2mat(dist_span);
if strcmp(opt.type,'max')
    target_dist=max(span_mat(all(span_mat>=0,2),:),[],2);
else
    target_dist=min(span_mat(all(span_mat>=0,2),:),[],2);
end
distci=bootci(1000,@(x) histcounts(x,0:7,'Normalization','probability'), target_dist)*100;
mm=histcounts(target_dist,0:7,'Normalization','probability')*100;
fh=figure('Color','w','Position',[100,100,240,180]);
hold on

bar(0:6,mm,'FaceColor','none','EdgeColor','k','LineWidth',1);
errorbar(0:6,mm,distci(1,:)-mm,distci(2,:)-mm,'k.','LineWidth',1);
ylabel('Normalized probability (%)')
set(gca,'XTick',0:6,'XTickLabel',7:-1:1,'FontSize',10);
if strcmp(opt.type,'max')
    xlabel(sprintf('Highest inter-struct. level, n=%d',numel(target_dist)));
else
    xlabel(sprintf('Lowest inter-struct. level, n=%d',numel(target_dist)));
end

ylim([0,100]);
exportgraphics(fh,fullfile('bzdata',sprintf('Ring_dist_prob_%s_%d.pdf',opt.type,msize)));
end
