function [metas,stats,fwd_rev,per_trial]=get_fc_coding(opt)
arguments
    opt.no_jitter (1,1) logical = true  % true for observed, false for jitter_shifted
    opt.shuffle (1,1) logical = false
    opt.fwd_rev (1,1) logical = false
    opt.per_trial (1,1) logical = false
end
persistent metas_ stats_ fwd_rev_ per_trial_ opt_
if isempty(metas_) || isempty(stats_) || isempty(fwd_rev_) || isempty(per_trial_) || ~isequaln(opt_,opt)
    sig=bz.load_sig_pair();
%     sess=unique(sig.sess);
    fl=dir(fullfile('fcdata','fc_coding_*.mat')); % data source jpsth\+bz\+fccoding\fc_coding_one_sess.m
    % {suids(fci,:),fwd_fc,fwd_fc-mean(fwd_shift,2),fwd_shift,rev_fc,rev_fc-mean(rev_shift,2),rev_shift}
    [metas,stats,fwd_rev]=deal([]);
    per_trial=cell(0);
    for fi=1:numel(fl)
        %     disp(fi);
        fstr=load(fullfile(fl(fi).folder,fl(fi).name));
        sess=fstr.onesess.fidx;
        % mapping FCSP to FC-meta
        sesssel=sig.sess==sess;
        ids=sig.suid(sesssel,:);
        mems=sig.mem_type(sesssel,:);
        regs=squeeze(sig.reg(sesssel,5,:));
        roots=squeeze(sig.reg(sesssel,1,:));
        %
        
        s1sel=fstr.onesess.trials(:,5)==4 & fstr.onesess.trials(:,8)==6 & all(fstr.onesess.trials(:,9:10),2);
        s2sel=fstr.onesess.trials(:,5)==8 & fstr.onesess.trials(:,8)==6 & all(fstr.onesess.trials(:,9:10),2);
        %
        if opt.shuffle
            pool=randsample([s1sel;s2sel],numel(s1sel)*2);
            s1sel=pool(1:numel(s1sel));
            s2sel=pool(numel(s1sel)+1:end);
        end
            
        e1sel=fstr.onesess.trials(:,5)==4 & fstr.onesess.trials(:,8)==6 & fstr.onesess.trials(:,10)==0;
        e2sel=fstr.onesess.trials(:,5)==8 & fstr.onesess.trials(:,8)==6 & fstr.onesess.trials(:,10)==0;
        if opt.no_jitter
            stat_idx=2;
        else
            stat_idx=3;
        end
        
        for fci=1:size(fstr.onesess.fc,1) %fc idx
            suid1=fstr.onesess.fc{fci,1}(1);
            suid2=fstr.onesess.fc{fci,1}(2);
            fcsel=ids(:,1)==suid1 & ids(:,2)==suid2;%map to meta
            mem=mems(fcsel,:);
            root=roots(fcsel,:);
            reg=regs(fcsel,:);
            if any(isempty(reg)) || any(reg==0,'all') || ~all(root==567,'all'), continue;end
            pertrl=fstr.onesess.fc{fci,stat_idx}>0;
            metas=[metas;sess,suid1,suid2,mem,reg];
            stats=[stats;mean(fstr.onesess.fc{fci,stat_idx}(s1sel)),mean(fstr.onesess.fc{fci,stat_idx}(s2sel)),...%1,2
                mean(fstr.onesess.fc{fci,stat_idx}(e1sel)),mean(fstr.onesess.fc{fci,stat_idx}(e2sel)),... %3,4
                nnz(pertrl(s1sel))/nnz(s1sel),nnz(pertrl(s2sel))/nnz(s2sel)];%5,6
            if opt.per_trial
                per_trial=[per_trial;{fstr.onesess.fc{fci,stat_idx}(s1sel),fstr.onesess.fc{fci,stat_idx}(s2sel),fstr.onesess.fc{fci,stat_idx}(e1sel),fstr.onesess.fc{fci,stat_idx}(e2sel)}];
            end
            if opt.fwd_rev
                fwd_rev=[fwd_rev;...
                    mean(fstr.onesess.fc{fci,2}(s1sel)>fstr.onesess.fc{fci,5}(s1sel)),...
                    mean(fstr.onesess.fc{fci,2}(s2sel)>fstr.onesess.fc{fci,5}(s2sel))];
            end
        end
        disp(size(metas,1));
    end
    metas_=metas;
    stats_=stats;
    fwd_rev_=fwd_rev;
    per_trial_=per_trial;
    opt_=opt;
else
    metas=metas_;
    stats=stats_;
    fwd_rev=fwd_rev_;
    per_trial=per_trial_;
end
