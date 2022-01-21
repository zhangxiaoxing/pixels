function fc_rate_svm()
arguments
    
end
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
load('OBM1Map.mat','OBM1map');
[metas,stats]=bz.fccoding.get_fc_coding_trial();
congru_diff=rpt_one('diff','congru',metas,stats,idmap,OBM1map);
keyboard()
%% could be very time consuming for large repeats!
%should backport from cluster.
parfor rpt=1:3
    [metas,stats,~]=bz.fccoding.get_fc_coding('no_jitter',true,'shuffle',true);
    out_shuf(rpt)=plot_one('shuffle',metas,stats,idmap,OBM1map,false);
end

exportgraphics(fh,'FC_coding_hier.pdf')
keyboard()
end

function out=rpt_one(rtype,mtype,meta_in,stat_in,idmap,index_map,opt)
arguments
    rtype (1,:) char {mustBeMember(rtype,{'same','diff','any'})}
    mtype (1,:) char {mustBeMember(mtype,{'congru','incong','nonmem'})}
    meta_in %sess,suid1,suid2,mem1,mem2,reg1,reg2
    stat_in
    idmap
    index_map
    opt.n_trial (1,1) double = 20
    opt.n_su (1,1) double = 100
end
% output -> 100x20x2 

w_key=arrayfun(@(x) index_map.isKey(char(idmap.ccfid2reg(x))),meta_in(:,6:7));
metas=meta_in(all(w_key,2),:); % only region in map kept below
stats=stat_in(all(w_key,2),:);%s1,s2,e1,e2

% index_value=arrayfun(@(x) index_map(char(idmap.ccfid2reg(x))),metas(:,6:7));
% for future separation of H2L and L2H

switch rtype
    case 'same'
        rsel = metas(:,6)==metas(:,7) ;
        
    case 'diff'
        congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2) & metas(:,6)~=metas(:,7) ;
        congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2) & metas(:,6)~=metas(:,7);
    case 'any' % only region in map kept, see above
        congrus1=metas(:,4)==2 & metas(:,5)==2 & all(~ismissing(stats),2);
        congrus2=metas(:,4)==4 & metas(:,5)==4 & all(~ismissing(stats),2);
end

switch mtype
    case 'congru'
        msel=all(ismember(metas(:,4:5),1:2),2) | all(ismember(metas(:,4:5),3:4),2);
    case 'incong'
        msel=any(ismember(metas(:,4:5),1:2),2) & any(ismember(metas(:,4:5),3:4),2);
    case 'nonmem' % only region in map kept, see above
        msel=all(metas(:,4:5)==0,2);
end




if ~strcmp(rtype,'shuffle')
    dists1=arrayfun(@(x) -diff(stats(x,1:2),1,2)./sum(stats(x,1:2),2),find(congrus1 & sum(stats(:,1:2),2)>10));
    dists2=arrayfun(@(x) diff(stats(x,1:2),1,2)./sum(stats(x,1:2),2),find(congrus2 & sum(stats(:,1:2),2)>10));
    dists1e=arrayfun(@(x) -diff(stats(x,3:4),1,2)./sum(stats(x,3:4),2),find(congrus1 & sum(stats(:,3:4),2)>10));
    dists2e=arrayfun(@(x) diff(stats(x,3:4),1,2)./sum(stats(x,3:4),2),find(congrus2& sum(stats(:,3:4),2)>10));
    
    mm=[mean([dists1;dists2]),mean([dists1e;dists2e])];
    cic=bootci(1000,@(x) mean(x),[dists1;dists2]);
    cie=bootci(1000,@(x) mean(x),[dists1e;dists2e]);
    out.mm=mm;
    out.cic=cic;
    out.cie=cie;
    out.wrsp=ranksum([dists1;dists2],[dists1e;dists2e]);
    out.datac=[dists1;dists2];
    out.datae=[dists1e;dists2e];
else
    dists1=arrayfun(@(x) -diff(stats(x,1:2),1,2)./sum(stats(x,1:2),2),find(congrus1 & sum(stats(:,1:2),2)>10));
    dists2=arrayfun(@(x) diff(stats(x,1:2),1,2)./sum(stats(x,1:2),2),find(congrus2 & sum(stats(:,1:2),2)>10));
    mm=mean([dists1;dists2]);
    out.mm=mm;
    out.datac=[dists1;dists2];
end

end