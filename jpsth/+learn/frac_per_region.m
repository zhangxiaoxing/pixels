if false
global_init;
%%%%%%%%%% fraction
wt_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","WT","load_file",false,"skip_stats",true);
wt_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','WT','extend6s',true);
ln_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Learning","load_file",false,"skip_stats",true);
ln_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);

nv_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Naive","load_file",false,"skip_stats",true);
nv_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Naive','extend6s',true);

[wt_map,wt_fh]=ephys.pct_reg_bars(wt_su_meta,wt_sel_meta,'xyscale',{'linear','linear'},'only_odor',true,'criteria','WT','skip_plot',true); % only need map_cells for tcom-frac corr
[ln_map,ln_fh]=ephys.pct_reg_bars(ln_su_meta,ln_sel_meta,'xyscale',{'linear','linear'},'only_odor',true,'criteria','Learning','skip_plot',true); % only need map_cells for tcom-frac corr
[nv_map,nv_fh]=ephys.pct_reg_bars(nv_su_meta,nv_sel_meta,'xyscale',{'linear','linear'},'only_odor',true,'criteria','Naive','skip_plot',false); % only need map_cells for tcom-frac corr
end
[sfrac,sidx]=sortrows(cell2mat(wt_map.olf.values.'),[-1,3]);
regs=subsref(wt_map.olf.keys,substruct('()',{sidx}));
frac_mat=[sfrac,nan(size(sfrac)+[0,3])];
chisqp=nan(size(frac_mat,1),1);
for lnkey=reshape(ln_map.olf.keys,1,[])
    [in,pos]=ismember(lnkey,regs);
    if in
        frac_mat(pos,4:6)=ln_map.olf(lnkey{1});
        % [~,p]=fishertest([frac_mat(pos,2),frac_mat(pos,3)-frac_mat(pos,2);frac_mat(pos,5),frac_mat(pos,6)-frac_mat(pos,5)]);
        [~,~,p]=crosstab([zeros(frac_mat(pos,3),1);ones(frac_mat(pos,6),1)],[(1:frac_mat(pos,3))>frac_mat(pos,2),(1:frac_mat(pos,6))>frac_mat(pos,5)]);
        chisqp(pos)=p;
    end
end

for nvkey=reshape(nv_map.olf.keys,1,[])
    [in,pos]=ismember(nvkey,regs);
    if in
        frac_mat(pos,7:9)=nv_map.olf(nvkey{1});
    end
end



return
fh=figure('Position',[100,100,1280,480]);
finisel=find(isfinite(frac_mat(:,4)));
bh=bar(frac_mat(finisel,[7,4,1]).*100,'grouped');
set(gca,'XTick',1:numel(finisel),'XTickLabel',regs(finisel));
% text(find(isnan(frac_mat(:,4))),repmat(10,nnz(isnan(frac_mat(:,4))),1),repmat('/',nnz(isnan(frac_mat(:,4))),1),'HorizontalAlignment','right')
for pp=1:numel(finisel)
    if chisqp(finisel(pp))<0.001
        text(pp,85,'***','HorizontalAlignment','center')
    elseif chisqp(finisel(pp))>0.05
        text(pp,85,'ns','HorizontalAlignment','center')
    else
        text(pp,85,sprintf('%.4f',chisqp(finisel(pp))),'HorizontalAlignment','center')
    end
end

ylabel('Percentage of memory neuron %')
ylim([0,100]);
legend(bh,{'Naive','Learning','Well-trained'});
savefig(fh,fullfile('binary','LN_memory_su_region_dist.fig'))

