function reg_avg_com_map=reg_avg_tcom_map(pct_meta,com_map,opt)

arguments
    pct_meta
    com_map
    opt.single_mod_thresh (1,1) double = 4
end

% TODO: proper encapsulation 
% eff_meta=ephys.effect_size_meta();
% sens_efsz=max(abs(eff_meta.cohen_d_olf),[],2);
% sens_win=[min(sens_efsz)./2,prctile(sens_efsz,[20:20:100])];
% 
% dur_efsz=max(abs(eff_meta.cohen_d_dur),[],2);
% dur_win=[min(dur_efsz)./2,prctile(dur_efsz,[20:20:100])];
% pct_meta=pct.get_pct_meta(eff_meta,sens_efsz,sens_win,dur_efsz,dur_win,'single_mod_thresh',opt.single_mod_thresh);
% com_map=wave.get_pct_com_map(pct_meta,'curve',true);
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);

fns=fieldnames(com_map);

stats=cell(0);

for fn=reshape(fns,1,[])
    sessid=str2double(replace(fn{1},'s',''));
    prefs=fieldnames(com_map.(fn{1}));
    for pp=reshape(prefs,1,[])
        suid=com_map.(fn{1}).(pp{1}).com.keys();
        com=com_map.(fn{1}).(pp{1}).com.values();
        % sessid,suid,com,pref,reg
        currstat=nan(numel(suid),3);
        currstat(:,1)=sessid;
        currstat(:,2)=reshape(cell2mat(suid),[],1);
        currstat(:,3)=reshape(cell2mat(com),[],1);
        
        %reg
        sess_suid=meta.allcid(meta.sess==sessid);
        sess_reg=meta.reg_tree(5,meta.sess==sessid);
        [~,idx]=ismember(currstat(:,2),sess_suid);
        reg=reshape(sess_reg(idx),[],1);
        stats=[stats;num2cell(currstat),repmat(pp,numel(suid),1),reg];
    end
end

stats=stats(~ismissing(stats(:,5)),:);
grey_regs=ephys.getGreyRegs('range','grey');
ureg=intersect(unique(stats(:,5)),grey_regs);

reg_avg_com_map=containers.Map('KeyType','char','ValueType','any');

for rr=reshape(ureg,1,[])
    com_vec=cell2mat(stats(strcmp(rr,stats(:,5)),3));
    reg_avg_com_map(char(rr))=mean(com_vec);
end

end