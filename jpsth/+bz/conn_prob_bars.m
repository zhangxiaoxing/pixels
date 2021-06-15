% input data from
% [sig,pair]=bz.load_sig_pair('pair',true)

function conn_prob_bars(sig,pair,opt)
arguments
    sig (1,1) struct
    pair (1,1) struct
    opt.dist (1,1) double = 5
end

addpath('k:\code\align\')
% [sig.reg_dist,pair.reg_dist]=get_conn_dist(sig,pair);
sess_cnt=max(sig.sess);
same_stats=struct();
same_stats.nm_nm=nan(sess_cnt,1);
same_stats.congr=nan(sess_cnt,1);
same_stats.incon=nan(sess_cnt,1);
same_stats.mem_nm=nan(sess_cnt,1);
same_stats.nm_mem=nan(sess_cnt,1);
diff_stats=same_stats;
mm=cell(7,1);
ci=cell(7,1);
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
%same + diff
% for distbin=1:7
%     mm{distbin}=cellfun(@(x) nanmean(x)*100,cell5);
%     ci{distbin}=cell2mat(cellfun(@(x) bootci(1000,@(y) nanmean(y)*100,x),cell5,'UniformOutput',false));
% end
for fld=["nm_nm","congr","incon","mem_nm","nm_mem"]
    if isfield(onesess,fld)
        same_stats.sums.(fld).mm=nanmean(same_stats.(fld)).*100;
        same_stats.sums.(fld).ci=bootci(1000,@(x) nanmean(x),same_stats.(fld)).*100;
        diff_stats.sums.(fld).mm=nanmean(diff_stats.(fld)).*100;
        diff_stats.sums.(fld).ci=bootci(1000,@(x) nanmean(x),diff_stats.(fld)).*100;
    end
end

% finisel=all(cell2mat(cellfun(@(x) ...
%     all(isfinite(same_stats.(x)),2)...
%     & all(isfinite(diff_stats.(x)),2),...
%     {'nm_nm','congr','incon','mem_nm','nm_mem'},'UniformOutput',false)),2);
% p=anova2([same_stats.nm_nm(finisel,:);same_stats.congr(finisel,:)],nnz(finisel),'off');
% disp(p);
colors={'k','c','m','b','r'};
fh=figure('Color','w','Position',[100,100,250,250]);
subplot(1,2,1);
hold on;
xidx=1;
for fld=["nm_nm","congr","incon","mem_nm","nm_mem"]
    if isfield(onesess,fld)
        bar(xidx,same_stats.sums.(fld).mm,'FaceColor','w','EdgeColor','k');
        errorbar(xidx,...
            same_stats.sums.(fld).mm,...
            same_stats.sums.(fld).ci(1)-same_stats.sums.(fld).mm,...
            same_stats.sums.(fld).ci(2)-same_stats.sums.(fld).mm,...
            '-k.');
        xidx=xidx+1;
    end
end
ylabel('Coupling fraction (%)');
set(gca(),'XTick',1:5,'XTickLabel',["Non-Non","Congru","Incongru","Mem-Non","Non-Mem"],...
    'XTickLabelRotation',60,'TickLabelInterpreter','none','FontSize',10,...
    'YTick',0:2)
subplot(1,2,2);
hold on;
xidx=1;
for fld=["nm_nm","congr","incon","mem_nm","nm_mem"]
    if isfield(onesess,fld)
        bar(xidx,diff_stats.sums.(fld).mm,'FaceColor','w','EdgeColor','k');
        errorbar(xidx,...
            diff_stats.sums.(fld).mm,...
            diff_stats.sums.(fld).ci(1)-diff_stats.sums.(fld).mm,...
            diff_stats.sums.(fld).ci(2)-diff_stats.sums.(fld).mm,...
            '-k.');
        xidx=xidx+1;
    end
end
ylim([0,1]);
set(gca(),'XTick',1:5,'XTickLabel',["Non-Non","Congru","Incongru","Mem-Non","Non-Mem"],...
    'XTickLabelRotation',60,'TickLabelInterpreter','none','FontSize',10,...
    'YTick',0:0.5:1)

% exportgraphics
% 
% arrayfun(@(x) errorbar(x,mm(x),ci{x}(1)-mm(x),ci{x}(2)-mm(x),'k.','LineWidth',1),1:5);
% set(gca,'XTick',1:5,'XTickLabel',{'NonMem-NonMem','NonMem-Mem','Mem-NonMem','Incongruent','Congruent'},'XTickLabelRotation',45);
% ylabel('Connection fraction (%)')
% xlim([0.5,5.5]);
% % exportgraphics(fh,fullfile('bzdata','conn_frac.pdf'));
end