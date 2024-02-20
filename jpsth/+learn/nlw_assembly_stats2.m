
% wt/ln/nv per session -> sampling su -> build essemble
% sample 100 for session > 200 su
% TODO: investigate sessions with su but no SC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global_init;
wt_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","WT","load_file",false,"skip_stats",true);
nv_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Naive","load_file",false,"skip_stats",true);
ln_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Learning","load_file",false,"skip_stats",true);

wt_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','WT');
nv_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','Naive');
ln_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','Learning');

wt_sig=bz.load_sig_sums_conn_file('pair',false,'inhibit',false,'criteria','WT');
ln_sig=bz.load_sig_sums_conn_file('pair',false,'inhibit',false,'criteria','Learning');
nv_sig=bz.load_sig_sums_conn_file('pair',false,'inhibit',false,'criteria','Naive');


SAMP_SU=100;
MIN_SESS_CNT=200;
PER_SESS_RPT=100;


critset={'WT','Learning','Naive'};
outtbl=table('Size',[0,4],'VariableTypes',{'categorical','int32','int32','cell'});


usess=intersect(unique(wt_su_meta.sess), unique(wt_sig.sess));

for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("WT "+num2str(sidx))
    end
    sess_sel=wt_su_meta.sess==sidx;
    if false
        s1sel=sess_sel & wt_sel_meta.wave_id==5;
        s2sel=sess_sel & wt_sel_meta.wave_id==6;
        nmsel=sess_sel & wt_sel_meta.wave_id==0;

        s1lbl=arrayfun(@(x) num2str(x),wt_su_meta.allcid(s1sel),'UniformOutput',false);
        s2lbl=arrayfun(@(x) num2str(x),wt_su_meta.allcid(s2sel),'UniformOutput',false);
        nmlbl=arrayfun(@(x) num2str(x),wt_su_meta.allcid(nmsel),'UniformOutput',false);
    end

    cnt=nnz(sess_sel);
    if cnt>MIN_SESS_CNT
        for rpt=1:PER_SESS_RPT
            samp_su=randsample(wt_su_meta.allcid(sess_sel),SAMP_SU,false);
            sig_sel=wt_sig.sess==sidx;
            sesssig=wt_sig.suid(sig_sel,:);
            samp_sig_sel=all(ismember(sesssig,samp_su),2);
            gh=digraph(string(sesssig(samp_sig_sel,1)),string(sesssig(samp_sig_sel,2)));
            [~,binsize]=gh.conncomp('Type','weak');
            % add non-mem congruent stats
            if false
                s1subg=gh.subgraph(intersect(gh.Nodes.Name,s1lbl));
                s2subg=gh.subgraph(intersect(gh.Nodes.Name,s2lbl));
                nmsubg=gh.subgraph(intersect(gh.Nodes.Name,nmlbl));
            end

            outtbl=[outtbl;cell2table({categorical({'WT'},critset),sidx,rpt,{binsize}})];
        end
    end
end


usess=intersect(unique(ln_su_meta.sess), unique(ln_sig.sess));
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Learning "+num2str(sidx))
    end
    sess_sel=ln_su_meta.sess==sidx;
    cnt=nnz(sess_sel);
    if cnt>MIN_SESS_CNT
        for rpt=1:PER_SESS_RPT
            samp_su=randsample(ln_su_meta.allcid(sess_sel),SAMP_SU,false);
            sig_sel=ln_sig.sess==sidx;
            sesssig=ln_sig.suid(sig_sel,:);
            samp_sig_sel=all(ismember(sesssig,samp_su),2);
            gh=digraph(string(sesssig(samp_sig_sel,1)),string(sesssig(samp_sig_sel,2)));
            [~,binsize]=gh.conncomp('Type','weak');
            outtbl=[outtbl;cell2table({categorical({'Learning'},critset),sidx,rpt,{binsize}})];
        end
    end
end

usess=intersect(unique(nv_su_meta.sess), unique(nv_sig.sess));
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Naive "+num2str(sidx))
    end
    sess_sel=nv_su_meta.sess==sidx;
    cnt=nnz(sess_sel);
    if cnt>MIN_SESS_CNT
        for rpt=1:PER_SESS_RPT
            samp_su=randsample(nv_su_meta.allcid(sess_sel),SAMP_SU,false);
            sig_sel=nv_sig.sess==sidx;
            sesssig=nv_sig.suid(sig_sel,:);
            samp_sig_sel=all(ismember(sesssig,samp_su),2);
            gh=digraph(string(sesssig(samp_sig_sel,1)),string(sesssig(samp_sig_sel,2)));
            [~,binsize]=gh.conncomp('Type','weak');
            
            outtbl=[outtbl;cell2table({categorical({'Naive'},critset),sidx,rpt,{binsize}})];
        end
    end
end
outtbl.Properties.VariableNames={'Criteria','Session','rpt','Size'};

% per-component-size
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=outtbl.Size{ii};
    kwv=[kwv;ss(ss>2).'];
    kwg=[kwg;repmat(outtbl.Criteria(ii),nnz(ss>2),1)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Connected component size');
title("p = "+num2str(p,4))


% largest-component-size
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=outtbl.Size{ii};
    kwv=[kwv;max(ss)];
    kwg=[kwg;outtbl.Criteria(ii)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Largest connected-component size');
title("p = "+num2str(p,4))
% anovan(kwv,{kwg})


% neurons involved in 3+ networks
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=outtbl.Size{ii};
    kwv=[kwv;sum(ss(ss>2))];
    kwg=[kwg;outtbl.Criteria(ii)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Percentage of neurons in ensembles of 3+ nodes (%)');
title("p = "+num2str(p,4))

% fragmentation

kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=outtbl.Size{ii};
    kwv=[kwv;nnz(ss>2)];
    kwg=[kwg;outtbl.Criteria(ii)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Unconnected components (fragments)');
title("p = "+num2str(p,4))



