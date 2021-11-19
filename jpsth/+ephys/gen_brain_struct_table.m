idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
meta=ephys.util.load_meta();
ctxsel=strcmp(meta.reg_tree(2,:),'CTX');
graysel=strcmp(meta.reg_tree(1,:),'CH') | strcmp(meta.reg_tree(1,:),'BS');
uregctx7=unique(meta.reg_tree(5,ctxsel));
ureggrey7=unique(meta.reg_tree(5,graysel & ~ctxsel));
%% table 1 brain structure in dataset
currid=1;
% tbl={'Index#','Abbreviation','Name','AllenCCF id','Structure tree path'};
tbl={'Index#','Abbreviation','Name','AllenCCF id','Single unit count'};
for rr=reshape([uregctx7,ureggrey7],1,[])
    if isempty(rr{1}), continue;end

    % with structural tree path
    %     tbl(end+1,:)={num2str(currid),rr{1},char(idmap.reg2full(rr{1})),num2str(idmap.reg2ccfid(rr{1})),...
    %         replace(strjoin(idmap.reg2tree(rr{1}),'-'),'root-grey-','')};

    tbl(end+1,:)={num2str(currid),rr{1},char(idmap.reg2full(rr{1})),num2str(idmap.reg2ccfid(rr{1})),...
        num2str(nnz(strcmp(meta.reg_tree(5,:),rr{1})))};

    currid=currid+1;
end
writecell(tbl,'Brain_structure.csv');

%% table 2 for other level of 

ureg1_6=unique([meta.reg_tree(1,graysel),meta.reg_tree(2,graysel),meta.reg_tree(3,graysel),meta.reg_tree(4,graysel)]);
currid=1;
tbl={'Index#','Abbreviation','Name','AllenCCF id','Structure tree path'};
for rr=reshape(ureg1_6,1,[])
    if isempty(rr{1}), continue;end
    tbl(end+1,:)={num2str(currid),rr{1},char(idmap.reg2full(rr{1})),num2str(idmap.reg2ccfid(rr{1})),...
        replace(strjoin(idmap.reg2tree(rr{1}),'-'),'root-grey-','')};
    currid=currid+1;
end
writecell(tbl,'Brain_structure_upper_level.csv');
