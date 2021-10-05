function [metas,stats,fwd_rev]=get_fc_coding(opt)
arguments
    opt.keep_trial (1,1) logical = false
    opt.no_shift (1,1) logical = false
end
persistent metas_ stats_ fwd_rev_ no_shift
if isempty(metas_) || isempty(stats_) || isempty(fwd_rev_) || no_shift~=opt.no_shift
    sig=bz.load_sig_pair();
    sess=unique(sig.sess);
    fl=dir('fcdata\fc_coding_*.mat');
    % {suids(fci,:),fwd_fc,fwd_fc-mean(fwd_shift,2),fwd_shift,rev_fc,rev_fc-mean(rev_shift,2),rev_shift}
    metas=[];
    stats=[];
    fwd_rev=[];
    for fi=1:numel(fl)
        %     disp(fi);
        fstr=load(fullfile(fl(fi).folder,fl(fi).name));
        sess=fstr.onesess.fidx;
        sesssel=sig.sess==sess;
        ids=sig.suid(sesssel,:);
        mems=sig.mem_type(sesssel,:);
        regs=squeeze(sig.reg(sesssel,5,:));
        roots=squeeze(sig.reg(sesssel,1,:));
        
        s1sel=fstr.onesess.trials(:,5)==4 & fstr.onesess.trials(:,8)==6 & all(fstr.onesess.trials(:,9:10),2);
        s2sel=fstr.onesess.trials(:,5)==8 & fstr.onesess.trials(:,8)==6 & all(fstr.onesess.trials(:,9:10),2);
        e1sel=fstr.onesess.trials(:,5)==4 & fstr.onesess.trials(:,8)==6 & fstr.onesess.trials(:,10)==0;
        e2sel=fstr.onesess.trials(:,5)==8 & fstr.onesess.trials(:,8)==6 & fstr.onesess.trials(:,10)==0;
        if opt.no_shift
            stat_idx=2;
        else
            stat_idx=3;
        end
        
        for fci=1:size(fstr.onesess.fc,1) %fc idx
            suid1=fstr.onesess.fc{fci,1}(1);
            suid2=fstr.onesess.fc{fci,1}(2);
            fcsel=ids(:,1)==suid1 & ids(:,2)==suid2;
            mem=mems(fcsel,:);
            root=roots(fcsel,:);
            reg=regs(fcsel,:);
            if any(isempty(reg)) || any(reg==0,'all') || ~all(root==567,'all'), continue;end
            pertrl=fstr.onesess.fc{fci,stat_idx}>0;
            metas=[metas;sess,suid1,suid2,mem,reg];
            stats=[stats;mean(fstr.onesess.fc{fci,stat_idx}(s1sel)),mean(fstr.onesess.fc{fci,stat_idx}(s2sel)),...%1,2
                mean(fstr.onesess.fc{fci,stat_idx}(e1sel)),mean(fstr.onesess.fc{fci,stat_idx}(e2sel)),... %3,4
                nnz(pertrl(s1sel))/nnz(s1sel),nnz(pertrl(s2sel))/nnz(s2sel)];
            fwd_rev=[fwd_rev;...
                mean(fstr.onesess.fc{fci,2}(s1sel)>fstr.onesess.fc{fci,5}(s1sel)),...
                mean(fstr.onesess.fc{fci,2}(s2sel)>fstr.onesess.fc{fci,5}(s2sel))];
        end
        disp(size(metas,1));
    end
    metas_=metas;
    stats_=stats;
    fwd_rev_=fwd_rev;
    no_shift=opt.no_shift;
else
    metas=metas_;
    stats=stats_;
    fwd_rev=fwd_rev_;
end
