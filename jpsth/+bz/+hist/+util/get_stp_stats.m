function [stats,memtypes]=get_stp_stats(ftick,prefix,opt)
arguments
    ftick (1,1) double {mustBeMember(ftick,[300,600,3000,6000])}
    prefix (1,:) char
    opt.suffix (1,:) char = []
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.any (1,1) logical = false
end
if ~isempty(opt.suffix) && ~startsWith(opt.suffix,'_'), opt.suffix=['_',opt.suffix];end
persistent stats_ ftick_ memtypes_ prefix_ suffix_ type_ criteria_

if isempty(stats_) || ftick~=ftick_ || ~strcmp(prefix,prefix_) || ~strcmp(opt.suffix,suffix_) || ~strcmp(opt.type,type_) || ~strcmp(criteria_,opt.criteria)
    disp('Updating STP stats from files');
    ftick_=ftick;
    fl=struct();
    
    if ~strcmp(opt.type,'MY')
        
        fl.congru=dir(fullfile('bzdata',sprintf('%s_stp_congru_*%d%s.mat',prefix,ftick,opt.suffix)));
        fl.incongru=dir(fullfile('bzdata',sprintf('%s_stp_incongru_*%d%s.mat',prefix,ftick,opt.suffix)));
        fl.nonmem=dir(fullfile('bzdata',sprintf('%s_stp_non-mem_*%d%s.mat',prefix,ftick,opt.suffix)));
    else
        if opt.any
            fl.any=dir(fullfile('mydata',sprintf('%s_stp_any_*%d%s.mat',prefix,ftick,opt.suffix)));
        else
            fl.congru=dir(fullfile('mydata',sprintf('%s_stp_congru_*%d%s.mat',prefix,ftick,opt.suffix)));
            fl.incongru=dir(fullfile('mydata',sprintf('%s_stp_incongru_*%d%s.mat',prefix,ftick,opt.suffix)));
            fl.nonmem=dir(fullfile('mydata',sprintf('%s_stp_non-mem_*%d%s.mat',prefix,ftick,opt.suffix)));
        end
    end
    
    memtypes=convertCharsToStrings(fieldnames(fl))';
    
    statfields=["postspk","skip","sess_suids"];
    stats=struct();
    
    for memtype=memtypes
        stats.(memtype)=struct();
        for sf=statfields, stats.(memtype).(sf)=[]; end
        stats.(memtype).sess=[];
        for fidx=1:size(fl.(memtype))
            fstr=load(fullfile(fl.(memtype)(fidx).folder,fl.(memtype)(fidx).name));
            %HOTFIX>>>>>>>>>>>
            if ~isfield(fstr,'skip')
                continue
            end
            if size(fstr.skip,2)>1
                fstr.skip=fstr.skip(1:size(fstr.sess_suids,1)).';
            end
            %<<<<<<<<<<<<<<<<<
            for sf=statfields, stats.(memtype).(sf)=[stats.(memtype).(sf);fstr.(sf)]; end
            stats.(memtype).sess=[stats.(memtype).sess;repmat(fstr.sess,size(fstr.sess_suids,1),1)];
        end
        stats.(memtype).reg=bz.hist.tag_hist_reg(stats.(memtype),'type',opt.type,'criteria',opt.criteria);
        [is_diff,is_same]=bz.hist.util.diff_at_level(stats.(memtype).reg);
        stats.(memtype).diff_reg=is_diff;
        stats.(memtype).same_reg=is_same;
    end
    stats_=stats;
    memtypes_=memtypes;
    prefix_=prefix;
    suffix_=opt.suffix;
    type_=opt.type;
    criteria_=opt.criteria;
else
    disp('Reusing STP stats in memory');
    stats=stats_;
    memtypes=memtypes_;
end

end