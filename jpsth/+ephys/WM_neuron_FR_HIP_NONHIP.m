load(fullfile('binary','su_meta.mat'));
load(fullfile('binary','wrs_mux_meta.mat'));
load(fullfile('binary','trials_dict.mat'),'trials_dict');

[gc,gr]=groupcounts(categorical(su_meta.reg_tree(5,ismember(wrs_mux_meta.wave_id,1:6)).'));
statreg=gr(gc>20);
statreg=cellstr(statreg(~isundefined(statreg)));

su_sums=cell2struct(cell(numel(statreg),1),statreg);
pref_samp=[4 4 8 8 4 8];
np_dur=[6 3 6 3 0 0];

for sess=reshape(unique(su_meta.sess),1,[])
    if ~any(ismember(su_meta.reg_tree(5,su_meta.sess==sess & ismember(wrs_mux_meta.wave_id,1:6)),statreg))
        continue
    end
    [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sess,'keep_trial',true,'jagged',true);
    trials=cell2mat(trials_dict(sess));

    for onereg=reshape(statreg,1,[])
        cids=su_meta.allcid(su_meta.sess==sess  & ismember(wrs_mux_meta.wave_id,1:6) & strcmp(su_meta.reg_tree(5,:),onereg).');
        for cid=reshape(cids,1,[])
            waveid=wrs_mux_meta.wave_id(su_meta.sess==sess & su_meta.allcid==cid & ismember(wrs_mux_meta.wave_id,1:6));
            pref_delay_trls=all(trials(:,9:10)==1,2) & trials(:,5)==pref_samp(waveid) & ismember(trials(:,8),setdiff([3 6],np_dur(waveid)));
            np_delay_trls=all(trials(:,9:10)==1,2) & trials(:,5)==setdiff([4,8],pref_samp(waveid)) & ismember(trials(:,8),setdiff([3 6],np_dur(waveid)));

            pref_sec=sum(trials(pref_delay_trls,8));
            np_sec=sum(trials(np_delay_trls,8));

        end
    end
end



