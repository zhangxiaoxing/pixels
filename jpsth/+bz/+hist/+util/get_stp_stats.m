function [stats,memtypes]=get_stp_stats(ftick,opt)
arguments
    ftick (1,1) double {mustBeMember(ftick,[300,600,3000,6000])}
    opt.prefix (1,:) char
    opt.suffix (1,:) char
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
end
if ~isempty(opt.suffix) && ~startsWith(opt.suffix,'_'), opt.suffix=['_',opt.suffix];end
persistent stats_ ftick_ memtypes_ prefix_ suffix_ type_


if isempty(stats_) || ftick~=ftick_ || ~strcmp(opt.prefix,prefix_) || ~strcmp(opt.suffix,suffix_) || ~strcmp(opt.type,type_)
    disp('Updating STP stats from files');
    ftick_=ftick;
    fl=struct();
    fl.congru=dir(fullfile('bzdata',sprintf('%s_stp_congru_*%d%s.mat',opt.prefix,ftick,opt.suffix)));
    fl.incongru=dir(fullfile('bzdata',sprintf('%s_stp_incongru_*%d%s.mat',opt.prefix,ftick,opt.suffix)));
    fl.nonmem=dir(fullfile('bzdata',sprintf('%s_stp_non-mem_*%d%s.mat',opt.prefix,ftick,opt.suffix)));
    memtypes=convertCharsToStrings(fieldnames(fl))';
    
    statfields=["fc_eff","fc_prob","postspk","skip","sess_suids"];
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
        stats.(memtype).reg=bz.hist.tag_hist_reg(stats.(memtype),'type',opt.type);
        [is_diff,is_same]=bz.hist.util.diff_at_level(stats.(memtype).reg);
        stats.(memtype).diff_reg=is_diff;
        stats.(memtype).same_reg=is_same;
    end
    stats_=stats;
    memtypes_=memtypes;
    prefix_=opt.prefix;
    suffix_=opt.suffix;
    type_=opt.type;
else
    disp('Reusing STP stats in memory');
    stats=stats_;
    memtypes=memtypes_;
end

end