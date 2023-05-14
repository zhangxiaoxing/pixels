% called in bz.fccoding.plot_coding

function [metas,stats]=get_fc_coding(sel_meta,opt)
arguments
    sel_meta
    opt.shuffle (1,1) logical = false
    opt.type (1,:) char {mustBeMember(opt.type,{'olf','dur'})} = 'olf'
    opt.correct_trials (1,1) double = 20
    opt.error_trials (1,1) double = 2
    opt.force_update (1,1) logical = true
    opt.incong (1,1) logical = false
end
persistent metas_ stats_ opt_

if isempty(metas_) || isempty(stats_) || ~isequaln(opt_,opt) || opt.force_update
    stat_idx=2; % 3 if deduce jitter
    sig=bz.load_sig_sums_conn_file();
    sig=bz.join_fc_waveid(sig,sel_meta.wave_id);

    fl=dir(fullfile('fccoding','fc_coding_*.mat')); % data source jpsth\+bz\+fccoding\fc_coding_one_sess.m
    % {suids(fci,:),fwd_fc,fwd_fc-mean(fwd_shift,2),fwd_shift,rev_fc,rev_fc-mean(rev_shift,2),rev_shift}
    stats=struct();
    [metas,stats.lbl1,stats.lbl2,stats.e1,stats.e2]=deal([]);
    for fi=1:numel(fl)

        fstr=load(fullfile(fl(fi).folder,fl(fi).name));
        sess=fstr.onesess.fidx;
        sesssel=sig.sess==sess;
        waveids=sig.waveid(sesssel,:);

        regs=squeeze(sig.reg(sesssel,5,:));
        %         roots=squeeze(sig.reg(sesssel,1,:));
        %
        switch opt.type
            case 'olf'
                if opt.incong
                    fcIdces=find(pct.su_pairs.get_incongru(waveids)); 
                else
                    congrusel=pct.su_pairs.get_congru(waveids); 
                    fcIdces=find(congrusel & all(ismember(waveids,1:6),2));
                end
                lbl1sel=fstr.onesess.trials(:,5)==4 & ismember(fstr.onesess.trials(:,8),[3 6]) & all(fstr.onesess.trials(:,9:10),2);
                lbl2sel=fstr.onesess.trials(:,5)==8 & ismember(fstr.onesess.trials(:,8),[3 6]) & all(fstr.onesess.trials(:,9:10),2);
                e1sel=fstr.onesess.trials(:,5)==4 & ismember(fstr.onesess.trials(:,8),[3 6]) & fstr.onesess.trials(:,10)==0;
                e2sel=fstr.onesess.trials(:,5)==8 & ismember(fstr.onesess.trials(:,8),[3 6]) & fstr.onesess.trials(:,10)==0;
            case 'dur'
                if opt.incong
                    fcIdces=find(pct.su_pairs.get_incongru(waveids));
                else
                    congrusel=pct.su_pairs.get_congru(waveids);
                    fcIdces=find(congrusel & all(ismember(waveids,[1:4,7:8]),2));
                end
                lbl1sel=ismember(fstr.onesess.trials(:,5),[4 8]) & fstr.onesess.trials(:,8)==3 & all(fstr.onesess.trials(:,9:10),2);
                lbl2sel=ismember(fstr.onesess.trials(:,5),[4 8]) & fstr.onesess.trials(:,8)==6 & all(fstr.onesess.trials(:,9:10),2);
                e1sel=ismember(fstr.onesess.trials(:,5),[4 8]) & fstr.onesess.trials(:,8)==3 & fstr.onesess.trials(:,10)==0;
                e2sel=ismember(fstr.onesess.trials(:,5),[4 8]) & fstr.onesess.trials(:,8)==6 & fstr.onesess.trials(:,10)==0;
        end
        

        if any([nnz(lbl1sel),nnz(lbl2sel)]<opt.correct_trials,'all')...
                || any([nnz(e1sel),nnz(e2sel)]<opt.error_trials,'all')
            continue
        end

        lbl1_indices=randsample(find(lbl1sel),opt.correct_trials);
        lbl2_indices=randsample(find(lbl2sel),opt.correct_trials);
        e1_indices=randsample(find(e1sel),opt.error_trials);
        e2_indices=randsample(find(e2sel),opt.error_trials);

        if opt.shuffle
            pool=randsample([lbl1sel;lbl2sel],numel(lbl1sel)*2);
            lbl1sel=pool(1:numel(lbl1sel));
            lbl2sel=pool(numel(lbl1sel)+1:end);
        end
            
        for fci=reshape(fcIdces,1,[]) %fc idx
            suids=reshape(fstr.onesess.fc{fci,1},1,[]);
            mem=waveids(fci,:); % TODO: waveid
%             root=roots(fcsel,:);
            reg=regs(fci,:);
            %TODO member of grey-regions, if needed
%             if any(isempty(reg)) || any(reg==0,'all') || ~all(root==567,'all'), continue;end
%             pertrl=fstr.onesess.fc{fci,stat_idx}>0;
            metas=[metas;sess,suids,mem,reg];
            stats.lbl1=[stats.lbl1,fstr.onesess.fc{fci,stat_idx}(lbl1_indices)];
            stats.lbl2=[stats.lbl2,fstr.onesess.fc{fci,stat_idx}(lbl2_indices)];
            stats.e1=[stats.e1,fstr.onesess.fc{fci,stat_idx}(e1_indices)];
            stats.e2=[stats.e2,fstr.onesess.fc{fci,stat_idx}(e2_indices)];
        end
%         disp(size(metas,1));
    end
    metas_=metas;
    stats_=stats;
    opt_=opt;
else
    metas=metas_;
    stats=stats_;
end
