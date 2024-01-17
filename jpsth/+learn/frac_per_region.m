global_init;
%%%%%%%%%% fraction
wt_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","WT","load_file",false,"skip_stats",true);
wt_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','WT','extend6s',true);
ln_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Learning","load_file",false,"skip_stats",true);
ln_sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);

[wt_map,wt_fh]=ephys.pct_reg_bars(wt_su_meta,wt_sel_meta,'xyscale',{'linear','linear'},'only_odor',true,'criteria','WT','skip_plot',true); % only need map_cells for tcom-frac corr
[ln_map,ln_fh]=ephys.pct_reg_bars(ln_su_meta,ln_sel_meta,'xyscale',{'linear','linear'},'only_odor',true,'criteria','Learning','skip_plot',true); % only need map_cells for tcom-frac corr

[sfrac,sidx]=sortrows(cell2mat(wt_map.olf.values.'),[-1,3]);
regs=subsref(wt_map.olf.keys,substruct('()',{sidx}));
frac_mat=[sfrac,nan(size(sfrac)+[0,1])];
for lnkey=reshape(ln_map.olf.keys,1,[])
    [in,pos]=ismember(lnkey,regs);
    if in
        frac_mat(pos,4:6)=ln_map.olf(lnkey{1});
        % [~,p]=fishertest([frac_mat(pos,2),frac_mat(pos,3)-frac_mat(pos,2);frac_mat(pos,5),frac_mat(pos,6)-frac_mat(pos,5)]);
        [~,~,p]=crosstab([zeros(frac_mat(pos,3),1);ones(frac_mat(pos,6),1)],[(1:frac_mat(pos,3))>frac_mat(pos,2),(1:frac_mat(pos,6))>frac_mat(pos,5)]);
        frac_mat(pos,7)=p;
    end
end

fh=figure('Position',[100,100,1280,480]);
finisel=find(isfinite(frac_mat(:,4)));
bh=bar(frac_mat(finisel,[4,1]).*100,'grouped');
set(gca,'XTick',1:numel(finisel),'XTickLabel',regs(finisel));
% text(find(isnan(frac_mat(:,4))),repmat(10,nnz(isnan(frac_mat(:,4))),1),repmat('/',nnz(isnan(frac_mat(:,4))),1),'HorizontalAlignment','right')
for pp=1:numel(finisel)
    if frac_mat(finisel(pp),7)<0.001
        text(pp,85,'***','HorizontalAlignment','center')
    elseif frac_mat(finisel(pp),7)>0.05
        text(pp,85,'ns','HorizontalAlignment','center')
    else
        text(pp,85,sprintf('%.4f',frac_mat(finisel(pp),7)),'HorizontalAlignment','center')
    end
end

ylabel('Percentage of memory neuron %')
ylim([0,100]);
legend(bh,{'Learning','Well-trained'});
savefig(fh,fullfile('binary','LN_memory_su_region_dist.fig'))

