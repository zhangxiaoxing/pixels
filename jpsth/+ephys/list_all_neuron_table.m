su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
uregs=setdiff(unique(su_meta.reg_tree(5,:)),{''});
reg_su_count=cellfun(@(x) nnz(strcmp(x,su_meta.reg_tree(5,:))),uregs);
[srt_cnt,srt_idx]=sort(reg_su_count,'descend');

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
%% per region stats

% tbl={'Index#','Abbreviation','Full name','Single unit count'};
tbl=cell(numel(srt_idx),7);
for currid=1:numel(srt_idx)
    regid=srt_idx(currid);
    % with structural tree path
    %     tbl(end+1,:)={num2str(currid),rr{1},char(idmap.reg2full(rr{1})),num2str(idmap.reg2ccfid(rr{1})),...
    %         replace(strjoin(idmap.reg2tree(rr{1}),'-'),'root-grey-','')};
    %     tbl(end+1,:)={num2str(currid),rr{1},char(idmap.reg2full(rr{1})),num2str(idmap.reg2ccfid(rr{1})),...
    %         num2str(nnz(strcmp(su_meta.reg_tree(5,:),rr{1})))};
    tbl(currid,:)={num2str(currid),uregs{regid},char(idmap.reg2full(uregs{regid})),num2str(srt_cnt(currid)),...
        num2str(nnz(strcmp(uregs{regid},su_meta.reg_tree(5,:)).' & ismember(wrs_mux_meta.wave_id,5:6))),...
        num2str(nnz(strcmp(uregs{regid},su_meta.reg_tree(5,:)).' & ismember(wrs_mux_meta.wave_id,1:4))),...
        num2str(nnz(strcmp(uregs{regid},su_meta.reg_tree(5,:)).' & ismember(wrs_mux_meta.wave_id,7:8)))};
end


tbl=[{'Index#','Abbreviation','Name','All neurons','Odor only','Encode both','Duration only'};tbl];
writecell(tbl,'Table_S1.csv');