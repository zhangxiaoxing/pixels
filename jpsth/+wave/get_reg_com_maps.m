function reg_com_maps=get_reg_com_maps(sel_meta,opt)
arguments
    sel_meta
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end

global_init;
com_map=wave.get_pct_com_map(sel_meta,'odor_only',true,'criteria',opt.criteria);

tcom3_maps=struct();
tcom6_maps=struct();
grey_regs=ephys.getGreyRegs('range','grey','criteria',opt.criteria);

[fcom3.odor_only.collection,fcom3.odor_only.com_meta]=wave.per_region_COM(...
    com_map,'sel_type','odor_only','com_field','com3','criteria',opt.criteria);
ureg=intersect(grey_regs,fcom3.odor_only.collection(:,2));
[~,tcidx]=ismember(ureg,fcom3.odor_only.collection(:,2));
tcom3_maps.odor_only=containers.Map(...
    ureg,num2cell(cellfun(@(x) x/4, fcom3.odor_only.collection(tcidx,1))));

[fcom6.odor_only.collection,fcom6.odor_only.com_meta]=wave.per_region_COM(...
    com_map,'sel_type','odor_only','com_field','com6','criteria',opt.criteria);
ureg=intersect(grey_regs,fcom6.odor_only.collection(:,2));
[~,tcidx]=ismember(ureg,fcom6.odor_only.collection(:,2));
tcom6_maps.odor_only=containers.Map(...
    ureg,num2cell(cellfun(@(x) x/4, fcom6.odor_only.collection(tcidx,1))));
reg_com_maps=cell2struct({tcom3_maps;tcom6_maps},{'tcom3_maps','tcom6_maps'});
end