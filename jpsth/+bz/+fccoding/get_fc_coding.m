function [metas,stats]=get_fc_coding(opt)
arguments
    opt.keep_trial (1,1) logical = false
end
persistent metas_ stats_
if isempty(metas_) || isempty(stats_)
    sig=bz.load_sig_pair();
    sess=unique(sig.sess);
    fl=dir('fcdata\fc_coding_*.mat');
    
    metas=[];
    stats=[];
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
        
        for ci=1:size(fstr.onesess.fc,1) %fc idx
            suid1=fstr.onesess.fc{ci,1}(1);
            suid2=fstr.onesess.fc{ci,1}(2);
            fcsel=ids(:,1)==suid1 & ids(:,2)==suid2;
            mem=mems(fcsel,:);
            root=roots(fcsel,:);
            reg=regs(fcsel,:);
            if any(isempty(reg)) || any(reg==0,'all') || ~all(root==567,'all'), continue;end
            pertrl=fstr.onesess.fc{ci,2}>mean(fstr.onesess.fc{ci,4},2);
            metas=[metas;sess,suid1,suid2,mem,reg];
            stats=[stats;mean(fstr.onesess.fc{ci,3}(s1sel)),mean(fstr.onesess.fc{ci,3}(s2sel)),...
                mean(fstr.onesess.fc{ci,3}(e1sel)),mean(fstr.onesess.fc{ci,3}(e2sel)),nnz(pertrl(s1sel))/nnz(s1sel),nnz(pertrl(s2sel))/nnz(s2sel)];
        end
        disp(size(metas,1));
    end
    metas_=metas;
    stats_=stats;
else
    metas=metas_;
    stats=stats_;
end
