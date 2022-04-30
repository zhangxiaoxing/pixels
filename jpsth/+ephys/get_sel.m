function [dur_mix,dur_exclu,dur_waveid,sens_exclu]=get_sel(opt)
arguments
    opt.ranksum_stats (1,1) logical = false
end

% persistent dur_exclu_ dur_mix_ dur_waveid_ sens_exclu_ opt_



    % clean house anova stats '22APR25
    % unfinished
    if isempty(dur_exclu_) || isempty (dur_mix_) || isempty(dur_waveid_) || ~isequaln(opt,opt_)
        anovameta=ephys.selectivity_anova();
        dur_waveid_=zeros(size(anovameta.sess));
        dur_all=any(anovameta.anovap(:,[2 4 6 7])<0.05,2);
        sens_all=any(anovameta.anovap(:,[1 4 5 7])<0.05,2);
        %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
        if opt.ranksum_stats
            meta=ephys.util.load_meta();
            waveid=ephys.get_wave_id(meta.sess,meta.allcid);
            dur_exclu_=dur_all & waveid==0; % relative more, may code odor to a degree
            dur_mix_=ismember(waveid,1:4); % relative less, odor coding double-checked
            sens_exclu_=ismember(waveid,5:6);
        else
            dur_exclu_=dur_all & ~sens_all; % relative less, removed interaction on first 3s
            sens_exclu_=sens_all & ~dur_all;
            dur_mix_=dur_all & sens_all & ~dur_exclu_ & ~sens_exclu_; %relative more, more than either-coding sense-sel su

        end
        dur_waveid_(dur_exclu_ & anovameta.dur_selidx>0)=3;
        dur_waveid_(dur_exclu_ & anovameta.dur_selidx<0)=6;
        opt_=opt;
    end

    dur_exclu=dur_exclu_;
    dur_mix=dur_mix_;
    dur_waveid=dur_waveid_;
    sens_exclu=sens_exclu_;
end
