function [stats,memtypes]=get_stp_stats(ftick,prefix,opt)
arguments
    ftick (1,1) double {mustBeMember(ftick,[300,600,3000,6000])}
    prefix (1,:) char
    opt.suffix (1,:) char = []
    opt.sesstype (1,:) char {mustBeMember(opt.sesstype,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.trialtype (1,:) char {mustBeMember(opt.trialtype,{'correct','error','any'})} = 'any'
    opt.datatype (1,:) char {mustBeMember(opt.datatype,{'postspk_out','fc_out'})}='postspk_out'
end
if ~isempty(opt.suffix) && ~startsWith(opt.suffix,'_'), opt.suffix=['_',opt.suffix];end

persistent stats_ ftick_ memtypes_ prefix_ suffix_ sesstype_ trialtype_ criteria_ datatype_

if isempty(stats_)...
        || ftick~=ftick_ ...
        || ~strcmp(prefix,prefix_) ...
        || ~strcmp(opt.suffix,suffix_) ...
        || ~strcmp(opt.sesstype,sesstype_) ...
        || ~strcmp(criteria_,opt.criteria) ...
        || ~strcmp(opt.trialtype,trialtype_) ...
        || ~strcmp(opt.datatype,datatype_) 
    disp('Updating STP stats from files');
    ftick_=ftick;
    fl=struct();
    
    if ~strcmp(opt.sesstype,'MY')
        fl.(opt.trialtype)=dir(fullfile('bzdata',sprintf('%s_stp_%s_*%d%s.mat',prefix,opt.trialtype,ftick,opt.suffix)));
    else
        fl.(opt.trialtype)=dir(fullfile('mydata',sprintf('%s_stp_%s_*%d%s.mat',prefix,opt.trialtype,ftick,opt.suffix)));
    end
    
    memtypes=convertCharsToStrings(fieldnames(fl))';
    
    statfields=["coeff","skip","sess_suids","pp","rsq"];
    stats=struct();
    
    for memtype=memtypes
        stats.(memtype)=struct();
        for sf=statfields, stats.(memtype).(sf)=[]; end
        stats.(memtype).sess=[];
        for fidx=1:size(fl.(memtype))
            fstrall=load(fullfile(fl.(memtype)(fidx).folder,fl.(memtype)(fidx).name));
            fstr=fstrall.(opt.datatype);
            fstr.sess_suids=fstrall.sess_suids;
            fstr.sess=fstrall.sess;
            if size(fstr.skip,2)>1
                fstr.skip=fstr.skip(1:size(fstr.sess_suids,1)).';
            end
            %<<<<<<<<<<<<<<<<<
            for sf=statfields, stats.(memtype).(sf)=[stats.(memtype).(sf);fstr.(sf)]; end
            stats.(memtype).sess=[stats.(memtype).sess;repmat(fstr.sess,size(fstr.sess_suids,1),1)];
        end
        [stats.(memtype).reg,stats.(memtype).mem_type]=bz.hist.tag_hist_meta(stats.(memtype),'type',opt.sesstype,'criteria',opt.criteria);
        [is_diff,is_same]=bz.hist.util.diff_at_level(stats.(memtype).reg);
        stats.(memtype).diff_reg=is_diff;
        stats.(memtype).same_reg=is_same;
    end
    stats_=stats;
    memtypes_=memtypes;
    prefix_=prefix;
    suffix_=opt.suffix;
    sesstype_=opt.sesstype;
    trialtype_=opt.trialtype;
    criteria_=opt.criteria;
    datatype_=opt.datatype;
else
    disp('Reusing STP stats in memory');
    stats=stats_;
    memtypes=memtypes_;
end

end