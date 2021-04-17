% Plots the history effect/short term plasticity/time constant between
% pairs of neuronal functional couplings. 
% Require dataset from \jpsth\+bz\+hist\hist_coeff_mem_nonmem.m
fl=struct();
fl.congru=dir(fullfile('bzdata','0331_stp_congru_*.mat'));
fl.incongru=dir(fullfile('bzdata','0331_stp_incongru_*.mat'));
fl.nonmem=dir(fullfile('bzdata','0331_stp_non-mem_*.mat'));

memtypes=convertCharsToStrings(fieldnames(fl))';
statfields=["fc_eff","fc_prob","postspk","maxiter","sess_suids"];
stats=struct();

for memtype=memtypes
    stats.(memtype)=struct();
    for sf=statfields, stats.(memtype).(sf)=[]; end
    stats.(memtype).sess=[];
    for fidx=1:size(fl.(memtype))
        fstr=load(fullfile(fl.(memtype)(fidx).folder,fl.(memtype)(fidx).name));
        for sf=statfields, stats.(memtype).(sf)=[stats.(memtype).(sf);fstr.(sf)]; end
        stats.(memtype).sess=[stats.(memtype).sess;repmat(fstr.sess,size(fstr.sess_suids,1),1)];
    end
    stats.(memtype).reg=bz.hist.tag_hist_reg(stats.(memtype));
end

%level3 ->{CTXpl,STR,TH,HY}
%level5 ->{PIR,AI,ORB,ILA,etc}
% stats_all=stats;
for level=[3,5]
    for memtype=memtypes
        [is_diff,is_same]=bz.hist.util.diff_at_level(stats.(memtype).reg,level);
        stats.(memtype).diff_reg=is_diff;
        stats.(memtype).same_reg=is_same;
    end
    fhspk_diff=bz.hist.plot_hist(stats.congru.postspk(stats.congru.maxiter(:,1)==0 & stats.congru.diff_reg,:),...
        stats.incongru.postspk(stats.incongru.maxiter(:,1)==0 & stats.incongru.diff_reg,:),...
        stats.nonmem.postspk(stats.nonmem.maxiter(:,1)==0 & stats.nonmem.diff_reg,:),...
        'title',sprintf('Diff. region branch %d', level));
    
    fhspk_same=bz.hist.plot_hist(stats.congru.postspk(stats.congru.maxiter(:,1)==0 & stats.congru.same_reg,:),...
        stats.incongru.postspk(stats.incongru.maxiter(:,1)==0 & stats.incongru.same_reg,:),...
        stats.nonmem.postspk(stats.nonmem.maxiter(:,1)==0 & stats.nonmem.same_reg,:),...
        'title',sprintf('Same region branch %d', level));
    
    fhfceff_diff=bz.hist.plot_hist(stats.congru.fc_eff(stats.congru.maxiter(:,2)==0 & stats.congru.diff_reg,:),...
        stats.incongru.fc_eff(stats.incongru.maxiter(:,2)==0 & stats.incongru.same_reg,:),...
        stats.nonmem.fc_eff(stats.nonmem.maxiter(:,2)==0 & stats.nonmem.diff_reg,:),...
        'type','fc_eff','title',sprintf('Diff. region branch %d', level));
    
    fhfceff_same=bz.hist.plot_hist(stats.congru.fc_eff(stats.congru.maxiter(:,2)==0 & stats.congru.same_reg,:),...
        stats.incongru.fc_eff(stats.incongru.maxiter(:,2)==0 & stats.incongru.same_reg,:),...
        stats.nonmem.fc_eff(stats.nonmem.maxiter(:,2)==0 & stats.nonmem.same_reg,:),...
        'type','fc_eff','title',sprintf('Same region branch %d', level));
    
    exportgraphics(fhspk_diff,sprintf('SPK_diff_L%d.pdf',level));
    exportgraphics(fhspk_same,sprintf('SPK_same_L%d.pdf',level));
    exportgraphics(fhfceff_diff,sprintf('FCEff_diff_L%d.pdf',level));
    exportgraphics(fhfceff_same,sprintf('FCEff_same_L%d.pdf',level));
end





fhspk=bz.hist.plot_hist(stats.congru.postspk(stats.congru.maxiter(:,1)==0,:),...
    stats.incongru.postspk(stats.incongru.maxiter(:,1)==0,:),...
    stats.nonmem.postspk(stats.nonmem.maxiter(:,1)==0,:));

fhfceff=bz.hist.plot_hist(stats.congru.fc_eff(stats.congru.maxiter(:,2)==0,:),...
    stats.incongru.fc_eff(stats.incongru.maxiter(:,2)==0,:),...
    stats.nonmem.fc_eff(stats.nonmem.maxiter(:,2)==0,:),'type','fc_eff');

fhfcprob=bz.hist.plot_hist(stats.congru.fc_prob(stats.congru.maxiter(:,3)==0,:),...
    stats.incongru.fc_prob(stats.incongru.maxiter(:,3)==0,:),...
    stats.nonmem.fc_prob(stats.nonmem.maxiter(:,3)==0,:),'type','fc_prob');
