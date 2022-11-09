function [metas,stats]=get_fc_coding(sel_meta,opt)
arguments
    sel_meta
    opt.shuffle (1,1) logical = false
    opt.type (1,:) char {mustBeMember(opt.type,{'olf','dur','mix'})} = 'olf'
    opt.correct_trials (1,1) double = 20
    opt.error_trials (1,1) double = 2
end
persistent metas_ stats_ opt_

if isempty(metas_) || isempty(stats_) || ~isequaln(opt_,opt)
    stat_idx=2; % 3 if deduce jitter
    sig=bz.load_sig_sums_conn_file();
    sig=bz.join_fc_waveid(sig,sel_meta.wave_id);

    fl=dir(fullfile('fccoding','fc_coding_*.mat')); % data source jpsth\+bz\+fccoding\fc_coding_one_sess.m
    % {suids(fci,:),fwd_fc,fwd_fc-mean(fwd_shift,2),fwd_shift,rev_fc,rev_fc-mean(rev_shift,2),rev_shift}
    [metas,stats]=deal([]);
    for fi=1:numel(fl)
        %     disp(fi);
        fstr=load(fullfile(fl(fi).folder,fl(fi).name));
        sess=fstr.onesess.fidx;
        % mapping FCSP to FC-meta
        sesssel=sig.sess==sess;
        waveids=sig.waveid(sesssel,:);
        congrusel=pct.su_pairs.get_congru(waveids); %TODO incongruent
        
%         ids=sig.suid(sesssel,:);
%         regs=squeeze(sig.reg(sesssel,5,:));
%         roots=squeeze(sig.reg(sesssel,1,:));
        % 
        switch opt.type
            case 'olf'
                fcIdces=find(congrusel & all(ismember(waveids,1:6),2));
                lbl1sel=fstr.onesess.trials(:,5)==4 & ismember(fstr.onesess.trials(:,8),[3 6]) & all(fstr.onesess.trials(:,9:10),2);
                lbl2sel=fstr.onesess.trials(:,5)==8 & ismember(fstr.onesess.trials(:,8),[3 6]) & all(fstr.onesess.trials(:,9:10),2);
                e1sel=fstr.onesess.trials(:,5)==4 & ismember(fstr.onesess.trials(:,8),[3 6]) & fstr.onesess.trials(:,10)==0;
                e2sel=fstr.onesess.trials(:,5)==8 & ismember(fstr.onesess.trials(:,8),[3 6]) & fstr.onesess.trials(:,10)==0;
        end
        

        if any([nnz(lbl1sel),nnz(lbl2sel)]<opt.correct_trials,'all')...
                || any([nnz(e1sel),nnz(e2sel)]<opt.error_trials,'all')
            continue
        end

        if opt.shuffle
            pool=randsample([lbl1sel;lbl2sel],numel(lbl1sel)*2);
            lbl1sel=pool(1:numel(lbl1sel));
            lbl2sel=pool(numel(lbl1sel)+1:end);
        end
            
        for fci=reshape(fcIdces,1,[]) %fc idx
%             suid1=fstr.onesess.fc{fci,1}(1);
%             suid2=fstr.onesess.fc{fci,1}(2);
%             fcsel=ids(:,1)==suid1 & ids(:,2)==suid2;%map to meta
%             mem=waveids(fcsel,:); % TODO: waveid
%             root=roots(fcsel,:);
%             reg=regs(fcsel,:);
            %TODO member of grey-regions, if needed
%             if any(isempty(reg)) || any(reg==0,'all') || ~all(root==567,'all'), continue;end
            pertrl=fstr.onesess.fc{fci,stat_idx}>0;
            metas=[metas;sess,suid1,suid2,mem,reg];
            stats=[stats;mean(fstr.onesess.fc{fci,stat_idx}(lbl1sel)),mean(fstr.onesess.fc{fci,stat_idx}(lbl2sel)),...%1,2
                mean(fstr.onesess.fc{fci,stat_idx}(e1sel)),mean(fstr.onesess.fc{fci,stat_idx}(e2sel)),... %3,4
                nnz(pertrl(lbl1sel))/nnz(lbl1sel),nnz(pertrl(lbl2sel))/nnz(lbl2sel)];%5,6

            %TODO:per-trial, label, fcsp
        

        end
        disp(size(metas,1));
    end
    metas_=metas;
    stats_=stats;
    opt_=opt;
else
    metas=metas_;
    stats=stats_;
end
