% Plots the history effect/short term plasticity/time constant between
% pairs of neuronal functional couplings. 
% Require dataset from \jpsth\+bz\+hist\hist_coeff_mem_nonmem.m
function plot_file_data(ftick,opt)
arguments
    ftick (1,1) double {mustBeMember(ftick,[300,600,3000,6000])}=3000;
    opt.prefix (1,:) char ='0428';
    opt.suffix (1,:) char ='laserOff';
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
end
[stats,~]=bz.hist.util.get_stp_stats(ftick,'prefix',opt.prefix,'suffix',opt.suffix,'type',opt.type);
binsize=ftick./30;
%level3 ->{CTXpl,STR,TH,HY}
%level5 ->{PIR,AI,ORB,ILA,etc}
% stats_all=stats;
for level=5

    fhspk_diff=bz.hist.plot_hist(stats.congru.postspk(stats.congru.skip(:,1)==0 & stats.congru.diff_reg(:,level),:),...
        stats.incongru.postspk(stats.incongru.skip(:,1)==0 & stats.incongru.diff_reg(:,level),:),...
        stats.nonmem.postspk(stats.nonmem.skip(:,1)==0 & stats.nonmem.diff_reg(:,level),:),...
        'title',sprintf('Cross region (%d ms bin) %s', binsize,replace(opt.suffix,'O',' o')),...
        'binw',binsize);
    
    fhspk_same=bz.hist.plot_hist(stats.congru.postspk(stats.congru.skip(:,1)==0 & stats.congru.same_reg(:,level),:),...
        stats.incongru.postspk(stats.incongru.skip(:,1)==0 & stats.incongru.same_reg(:,level),:),...
        stats.nonmem.postspk(stats.nonmem.skip(:,1)==0 & stats.nonmem.same_reg(:,level),:),...
        'title',sprintf('Within region (%d ms bin) %s', binsize,replace(opt.suffix,'O',' o')),...
        'binw',binsize);

    exportgraphics(fhspk_diff,sprintf('STP_diff_%dms_%s_%s_%s.pdf',binsize,opt.type,opt.prefix,opt.suffix));
    exportgraphics(fhspk_same,sprintf('STP_same_%dms_%s_%s_%s.pdf',binsize,opt.type,opt.prefix,opt.suffix));

end


