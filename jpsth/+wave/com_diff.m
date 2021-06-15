function comdiff_stats=com_diff(opt)

arguments
    opt.level (1,1) double {mustBeInteger} = 5
    opt.diff (1,1) logical = true
    opt.subplot (1,3) double {mustBeInteger,mustBePositive} = [1,1,1]
    opt.to_plot (1,1) logical = false
    opt.bin_edge (1,:) double = -2000:200:2000
    opt.peak (1,1) logical = false
end
com_map=wave.get_com_map('peak',opt.peak);
sig=bz.load_sig_pair('type','neupix','prefix','BZWT','criteria','WT');
[is_diff,is_same]=bz.util.diff_at_level(sig.reg);
comdiff_stats=[];
%mem_type [2 2] then [4 4]
if opt.diff
    trans1sel=find(all(sig.mem_type==2,2) & is_diff(:,opt.level)).';
    trans2sel=find(all(sig.mem_type==4,2) & is_diff(:,opt.level)).';
else
    trans1sel=find(all(sig.mem_type==2,2) & is_same(:,opt.level)).';
    trans2sel=find(all(sig.mem_type==4,2) & is_same(:,opt.level)).';
end

for ii=trans1sel
    sess=sig.sess(ii);
    suid=sig.suid(ii,:);
    com=arrayfun(@(x) com_map.(['s',num2str(sess)]).s1(x),suid);
    comdiff_stats=[comdiff_stats;diff(com)*250];
end

for ii=trans2sel
    sess=sig.sess(ii);
    suid=sig.suid(ii,:);
    com=arrayfun(@(x) com_map.(['s',num2str(sess)]).s2(x),suid);
    comdiff_stats=[comdiff_stats;diff(com)*250];
end

if opt.to_plot
    if isempty(opt.subplot)
        figure('Color','w');
    else
        subplot(opt.subplot(1),opt.subplot(2),opt.subplot(3))
    end
    histogram(comdiff_stats,opt.bin_edge,'Normalization','probability');
    ylabel('Probability');
    xlabel('Center of mass lag (post - pre, ms)');
    ylim([0,0.3]);
    if opt.diff,prefix='Between';else,prefix='Within';end
    title(sprintf('%s level %d region',prefix,opt.level));
    [~,p]=ttest(comdiff_stats);
    text(max(xlim()),max(ylim()),sprintf('bias=%.1fms, p=%0.3f',nanmean(comdiff_stats),p),'HorizontalAlignment','right','VerticalAlignment','top','FontSize',12);
end
end

