function [olf_sw,dur_sw,mux_sw]=get_switched(wrs_mux_meta,opt)
arguments
    wrs_mux_meta
    opt.odor_only
end

if opt.odor_only
    otypesel=ismember(wrs_mux_meta.wave_id,1:6);
    omask=wrs_mux_meta.p_olf<0.05 | wrs_mux_meta.p_mux<0.05;
    omask(~otypesel,:)=0;
    omasked=wrs_mux_meta.o_pref_id.*(omask);
    % olf_switched:bool
    olf_sw3=any(omasked==5,2) & any(omasked==6,2);

    otypesel=ismember(wrs_mux_meta.wave_id,[2 4 5 6]);
    omask=wrs_mux_meta.p_olf6<0.05;
    omask(~otypesel,:)=0;
    omasked=wrs_mux_meta.o_pref_id6.*(omask);
    % olf_switched:bool
    olf_sw6=any(omasked==5,2) & any(omasked==6,2);
    
    olf_sw=olf_sw3 | olf_sw6;
    [dur_sw,mux_sw]=deal(false(size(olf_sw)));
else

    otypesel=(ismember(wrs_mux_meta.wave_id,5:6));
    omask=wrs_mux_meta.p_olf<0.05;
    omask(~otypesel,:)=0;
    omasked=wrs_mux_meta.o_pref_id.*(omask);
    % olf_switched:bool
    olf_sw=any(ismember(omasked,1:2),2) & any(ismember(omasked,3:4),2);


    dtypesel=(ismember(wrs_mux_meta.wave_id,7:8));
    dmask=wrs_mux_meta.p_dur<0.05;
    dmask(~dtypesel,:)=0;
    dmasked=wrs_mux_meta.d_pref_id.*(dmask);
    % dur_switched:bool
    dur_sw=any(ismember(dmasked,[1 3]),2) & any(ismember(dmasked,[2 4]),2);

    % >>>>>>>>>>> different type of mux switched >>>>>>>>>>>>>>>>>>>>>>>
    mtypesel=(ismember(wrs_mux_meta.wave_id,1:4));
    mmmask=(wrs_mux_meta.p_mux<0.05);
    mmmask(~mtypesel,:)=0;
    mmmasked=wrs_mux_meta.m_pref_id.*(mmmask);
    mm_switched_sel=arrayfun(@(x)numel(unique(mmmasked(x,:)))-any(mmmasked(x,:)==0),1:size(mmmasked,1))>1;

    omtypesel=(ismember(wrs_mux_meta.wave_id,1:4));
    ommask=wrs_mux_meta.p_olf<0.05;
    ommask(~omtypesel,:)=0;
    ommasked=wrs_mux_meta.o_pref_id.*(ommask);
    om_switched_sel=any(ismember(ommasked,1:2),2) & any(ismember(ommasked,3:4),2);

    dmtypesel=(ismember(wrs_mux_meta.wave_id,1:4));
    dmmask=wrs_mux_meta.p_dur<0.05;
    dmmask(~dmtypesel,:)=0;
    dmmasked=wrs_mux_meta.d_pref_id.*(dmmask);
    dm_switched_sel=any(ismember(dmmasked,[1 3]),2) & any(ismember(dmmasked,[2 4]),2);
    % mux_switched:bool
    mux_sw=(dm_switched_sel | mm_switched_sel.' | om_switched_sel);
end
end
