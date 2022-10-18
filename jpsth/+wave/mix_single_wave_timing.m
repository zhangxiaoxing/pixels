global_init;
meta=ephys.util.load_meta('skip_stats',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true);


% function com_map2meta_like()
pref2waveid=containers.Map({'s1d3','s1d6','s2d3','s2d6','olf_s1','olf_s2','dur_d3','dur_d6'},...
    num2cell(1:8));
com_meta=struct();
[com_meta.sess,com_meta.allcid,com_meta.wave_id,com_meta.com]=deal([]);
sess=fieldnames(com_map);
for si=1:numel(sess)
    prefs=fieldnames(com_map.(sess{si}));
    for pi=1:numel(prefs)
        cmap=com_map.(sess{si}).(prefs{pi}).com;
        keys=cell2mat(cmap.keys());
        values=cell2mat(cmap.values());
        com_meta.sess=[com_meta.sess;...
            repmat(str2double(replace(sess{si},'s','')),...
            numel(keys),1)];
        com_meta.allcid=[com_meta.allcid;...
            reshape(keys,[],1)];
        com_meta.com=[com_meta.com;...
            reshape(values,[],1)];
        com_meta.wave_id=[com_meta.wave_id;...
            repmat(pref2waveid(prefs{pi}),numel(keys),1)];
    end
end

% end
