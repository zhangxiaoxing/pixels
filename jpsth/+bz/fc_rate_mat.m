if false
pthresh=gather_config.FC_thresh;
global_init();
wrs_mux_meta=ephys.get_wrs_mux_meta();
[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
sig=bz.join_fc_waveid(sig,wrs_mux_meta.wave_id);
pair=bz.join_fc_waveid(pair,wrs_mux_meta.wave_id);
greys=ephys.getGreyRegs('range','grey');
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
greys_id=int32(cell2mat(idmap.reg2ccfid.values(greys)));
sig.congrusel=pct.su_pairs.get_congru(sig.waveid);
pair.congrusel=pct.su_pairs.get_congru(pair.waveid);
% sig.reg_sel=all(ismember(sig.reg(:,5,:),greys_id),3);
% pair.reg_sel=all(ismember(pair.reg(:,5,:),greys_id),3);

sig.mix_sel=sig.congrusel & any(ismember(sig.waveid,1:4),2);
pair.mix_sel=pair.congrusel & any(ismember(pair.waveid,1:4),2);

sig.olf_sel=sig.congrusel & any(ismember(sig.waveid,5:6),2) & ~sig.mix_sel;
pair.olf_sel=pair.congrusel & any(ismember(pair.waveid,5:6),2) & ~pair.mix_sel;

sig.dur_sel=sig.congrusel & any(ismember(sig.waveid,7:8),2) & ~sig.mix_sel;
pair.dur_sel=pair.congrusel & any(ismember(pair.waveid,7:8),2) & ~pair.mix_sel;

both_ratio=nan(numel(greys_id),numel(greys_id),3);
olf_ratio=nan(numel(greys_id),numel(greys_id),3);
dur_ratio=nan(numel(greys_id),numel(greys_id),3);

for leadidx=1:numel(greys_id)
    lead_reg=greys_id(leadidx);
    for followidx=1:numel(greys_id)
        follow_reg=greys_id(followidx);
        preg_sel=(pair.reg(:,5,1)==lead_reg & pair.reg(:,5,2)==follow_reg);
% both
        pairn=nnz(preg_sel & pair.mix_sel);
        if pairn>pthresh
            sreg_sel=(sig.reg(:,5,1)==lead_reg & sig.reg(:,5,2)==follow_reg);
            sign=nnz(sreg_sel & sig.mix_sel);
            [hhat,cci]=binofit(sign,pairn);
            both_ratio(leadidx,followidx,1)=hhat;
            both_ratio(leadidx,followidx,2:3)=cci;
        end
% olf
        pairn=nnz(preg_sel & pair.olf_sel);
        if pairn>pthresh
            sreg_sel=(sig.reg(:,5,1)==lead_reg & sig.reg(:,5,2)==follow_reg);
            sign=nnz(sreg_sel & sig.olf_sel);
            [hhat,cci]=binofit(sign,pairn);
            olf_ratio(leadidx,followidx,1)=hhat;
            olf_ratio(leadidx,followidx,2:3)=cci;
        end
% dur
        pairn=nnz(preg_sel & pair.dur_sel);
        if pairn>pthresh
            sreg_sel=(sig.reg(:,5,1)==lead_reg & sig.reg(:,5,2)==follow_reg);
            sign=nnz(sreg_sel & sig.dur_sel);
            [hhat,cci]=binofit(sign,pairn);
            dur_ratio(leadidx,followidx,1)=hhat;
            dur_ratio(leadidx,followidx,2:3)=cci;
        end
    end
end

blame=vcs.blame();
save('fc_rate_mat.mat','blame','olf_ratio','dur_ratio','both_ratio','greys_id')
else
    load('fc_rate_mat.mat','olf_ratio','dur_ratio','both_ratio','greys_id')
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
end

% within
both_within=arrayfun(@(x) both_ratio(x,x,1),1:numel(greys_id));
olf_within=arrayfun(@(x) olf_ratio(x,x,1),1:numel(greys_id));
dur_within=arrayfun(@(x) dur_ratio(x,x,1),1:numel(greys_id));

[olf_within,olf_within_idx]=sort(olf_within,'descend','MissingPlacement','last');
[dur_within,dur_within_idx]=sort(dur_within,'descend','MissingPlacement','last');
[both_within,both_within_idx]=sort(both_within,'descend','MissingPlacement','last');

plotn=5;
figure()
tiledlayout(2,3)
nexttile()
bar(olf_within(1:plotn))
set(gca(),'XTick',1:plotn,'XTickLabelRotation',90,'XTickLabel',...
    string(idmap.ccfid2reg.values(num2cell(greys_id(olf_within_idx(1:plotn))))),...
    'YTick',0:0.04:0.08,'YTickLabel',0:4:8);
ylim([0,0.085])
ylabel('FC rate (%)')
title({'Olfactory','within'})
nexttile()
bar(dur_within(1:plotn))
set(gca(),'XTick',1:plotn,'XTickLabelRotation',90,'XTickLabel',...
    string(idmap.ccfid2reg.values(num2cell(greys_id(dur_within_idx(1:plotn))))),...
    'YTick',0:0.04:0.08,'YTickLabel',0:4:8);
ylim([0,0.085])
ylabel('FC rate (%)')
title({'Duration','within'})

nexttile()
bar(both_within(1:plotn))
set(gca(),'XTick',1:plotn,'XTickLabelRotation',90,'XTickLabel',...
    string(idmap.ccfid2reg.values(num2cell(greys_id(both_within_idx(1:plotn))))),...
    'YTick',0:0.04:0.08,'YTickLabel',0:4:8);
ylim([0,0.085])
ylabel('FC rate (%)')
title({'Both','within'})

% cross
both_cross=both_ratio(:,:,1);
olf_cross=olf_ratio(:,:,1);
dur_cross=dur_ratio(:,:,1);
for ii=1:numel(greys_id)
    both_cross(ii,ii)=NaN;
    olf_cross(ii,ii)=NaN;
    dur_cross(ii,ii)=NaN;
end

[olf_cross,olf_cross_idx]=sort(olf_cross(:),'descend','MissingPlacement','last');
[dur_cross,dur_cross_idx]=sort(dur_cross(:),'descend','MissingPlacement','last');
[both_cross,both_cross_idx]=sort(both_cross(:),'descend','MissingPlacement','last');
% figure()
% tiledlayout(1,3)
nexttile()
[row,col]=ind2sub([numel(greys_id),numel(greys_id)],olf_cross_idx(1:plotn));
xtk=arrayfun(@(x) string(idmap.ccfid2reg(greys_id(row(x))))+"-"+string(idmap.ccfid2reg(greys_id(col(x)))),1:plotn);
bar(olf_cross(1:plotn))
ylim([0,0.06])
set(gca(),'XTick',1:plotn,'XTickLabelRotation',90,'XTickLabel',xtk,...
    'YTick',0:0.02:0.08,'YTickLabel',0:2:8);
title({'Olfactory','cross'})

nexttile()
[row,col]=ind2sub([numel(greys_id),numel(greys_id)],dur_cross_idx(1:plotn));
xtk=arrayfun(@(x) string(idmap.ccfid2reg(greys_id(row(x))))+"-"+string(idmap.ccfid2reg(greys_id(col(x)))),1:plotn);
bar(dur_cross(1:plotn))
ylim([0,0.06])
set(gca(),'XTick',1:plotn,'XTickLabelRotation',90,'XTickLabel',xtk,...
    'YTick',0:0.02:0.08,'YTickLabel',0:2:8);
title({'Duration','cross'})

nexttile()
[row,col]=ind2sub([numel(greys_id),numel(greys_id)],both_cross_idx(1:plotn));
xtk=arrayfun(@(x) string(idmap.ccfid2reg(greys_id(row(x))))+"-"+string(idmap.ccfid2reg(greys_id(col(x)))),1:plotn);
bar(both_cross(1:plotn))
ylim([0,0.06])
set(gca(),'XTick',1:plotn,'XTickLabelRotation',90,'XTickLabel',xtk,...
    'YTick',0:0.02:0.08,'YTickLabel',0:2:8);
title({'Both','cross'})

