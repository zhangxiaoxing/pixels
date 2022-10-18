global_init;
sumeta=ephys.util.load_meta('skip_stats',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true);

% function com_map2meta_like()
pref2waveid=containers.Map({'s1d3','s1d6','s2d3','s2d6','olf_s1','olf_s2','dur_d3','dur_d6'},...
    num2cell(1:8));
com_meta=struct();
[com_meta.sess,com_meta.allcid,com_meta.wave_id,com_meta.com]=deal([]);
com_meta.reg_tree=cell(0);
sess=fieldnames(com_map);
for si=1:numel(sess)
    sessid=str2double(replace(sess{si},'s',''));
    sesscid=sumeta.allcid(sumeta.sess==sessid);
    sessreg=sumeta.reg_tree(:,sumeta.sess==sessid);
    prefs=fieldnames(com_map.(sess{si}));
    for pi=1:numel(prefs)
        cmap=com_map.(sess{si}).(prefs{pi}).com;
        keys=cell2mat(cmap.keys());
        values=cell2mat(cmap.values());
        com_meta.sess=[com_meta.sess;...
            repmat(sessid,numel(keys),1)];
        com_meta.allcid=[com_meta.allcid;...
            reshape(keys,[],1)];
        com_meta.com=[com_meta.com;...
            reshape(values,[],1)];
        com_meta.wave_id=[com_meta.wave_id;...
            repmat(pref2waveid(prefs{pi}),numel(keys),1)];

        [~,sidx]=ismember(uint16(keys),sesscid);
        com_meta.reg_tree=[com_meta.reg_tree,sessreg(:,sidx)];
    end
end

% end
if false
mean(com_meta.com(ismember(com_meta.wave_id,1:4)))
mean(com_meta.com(ismember(com_meta.wave_id,5:6)))
mean(com_meta.com(ismember(com_meta.wave_id,7:8)))

ranksum(com_meta.com(ismember(com_meta.wave_id,1:4)),...
    com_meta.com(ismember(com_meta.wave_id,5:6)))
end


mixcdf=histcounts(com_meta.com(ismember(com_meta.wave_id,1:4)),0:0.25:12,'Normalization','cdf');
olfcdf=histcounts(com_meta.com(ismember(com_meta.wave_id,5:6)),0:0.25:12,'Normalization','cdf');
durcdf=histcounts(com_meta.com(ismember(com_meta.wave_id,7:8)),0:0.25:12,'Normalization','cdf');

figure()
hold on
plot(mixcdf,'-k')
plot(olfcdf,'-r')
plot(durcdf,'-b')
ylabel('Cumulated density (%)')
set(gca(),'XTick',0:16:48,'XTickLabel',0:1:3,'YLim',[0,1])


figure()
histogram(com_meta.com(ismember(com_meta.wave_id,1:4)))
