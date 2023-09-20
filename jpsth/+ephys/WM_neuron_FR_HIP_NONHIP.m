load(fullfile('binary','su_meta.mat'));
load(fullfile('binary','wrs_mux_meta.mat'));

[gc,gr]=groupcounts(categorical(su_meta.reg_tree(5,ismember(wrs_mux_meta.wave_id,1:6)).'));
statreg=gr(gc>20);
statreg=cellstr(statreg(~isundefined(statreg)));

for sess=reshape(unique(su_meta.sess),1,[])
    ismember(su_meta.reg_tree(5,su_meta.sess==sess),statreg);

end

