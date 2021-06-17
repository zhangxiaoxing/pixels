% input data from
% [sig,pair]=bz.load_sig_pair('pair',true)

function conn_prob_bars(sig,pair,opt)
arguments
    sig (1,1) struct
    pair (1,1) struct
    opt.dist (1,1) double = 5
end

addpath('k:\code\align\')
sess_cnt=max(sig.sess);
same_stats=struct();
[same_stats.nm_nm,same_stats.congr,same_stats.incon,same_stats.mem_nm,same_stats.nm_mem]...
    =deal(nan(sess_cnt,1));
diff_stats=same_stats;
[sig_diff,sig_same]=bz.util.diff_at_level(sig.reg);
[pair_diff,pair_same]=bz.util.diff_at_level(pair.reg);

for i=1:sess_cnt
    if rem(i,20)==0, disp(i);end
    %same
    sess_sig_type=sig.mem_type(sig.sess==i & sig_same(:,opt.dist),:);
    sess_pair_type=pair.mem_type(pair.sess==i & pair_same(:,opt.dist),:);
    onesess=bz.bars_util.get_ratio(sess_sig_type,sess_pair_type);
    for fld=["nm_nm","congr","incon","mem_nm","nm_mem"]
        if isfield(onesess,fld)
            same_stats.(fld)(i)=onesess.(fld);
        end
    end
    %diff
    sess_sig_type=sig.mem_type(sig.sess==i & sig_diff(:,opt.dist),:);
    sess_pair_type=pair.mem_type(pair.sess==i & pair_diff(:,opt.dist),:);
    onesess=bz.bars_util.get_ratio(sess_sig_type,sess_pair_type);
    for fld=["nm_nm","congr","incon","mem_nm","nm_mem"]
        if isfield(onesess,fld)
            diff_stats.(fld)(i)=onesess.(fld);
        end
    end
end
flds=["nm_nm","congr","incon","mem_nm","nm_mem"];
samemat=cell2mat(arrayfun(@(x) same_stats.(x),flds,'UniformOutput',false));
diffmat=cell2mat(arrayfun(@(x) diff_stats.(x),flds,'UniformOutput',false));
finisel=all(isfinite([samemat,diffmat]),2);

same_stats.sums.mm=mean(samemat(finisel,:)).*100;
same_stats.sums.ci=bootci(1000,@(x) mean(x),samemat(finisel,:)).*100;
diff_stats.sums.mm=mean(diffmat(finisel,:)).*100;
diff_stats.sums.ci=bootci(1000,@(x) mean(x),diffmat(finisel,:)).*100;

psame=anova1(samemat(finisel,:),flds,'off');
pdiff=anova1(diffmat(finisel,:),flds,'off');

fh=figure('Color','w','Position',[100,100,250,250]);
subplot(1,2,1);hold on;
plotOne(same_stats,psame)
subplot(1,2,2);hold on;
plotOne(diff_stats,pdiff)
end


function plotOne(data,p)
bar(1:5,data.sums.mm,'FaceColor','w','EdgeColor','k');
errorbar(1:5,...
    data.sums.mm,...
    data.sums.ci(1,:)-data.sums.mm,...
    data.sums.ci(2,:)-data.sums.mm,...
    'k.');

if max(ylim())<1
    ylim([0,1]);
end
    
set(gca(),'XTick',1:5,'XTickLabel',["Non-Non","Congru","Incongru","Mem-Non","Non-Mem"],...
    'XTickLabelRotation',60,'TickLabelInterpreter','none','FontSize',10,...
    'YTick',0:0.5:max(ylim()))
text(max(xlim()),max(ylim()),sprintf('p=%.3f',p),'HorizontalAlignment','right','VerticalAlignment','top');
end