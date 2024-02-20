%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% congru incongru nonmemory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% wt/ln/nv per session -> sampling su -> build essemble
% sample 15 for each selectivity type > 20 su
% TODO: investigate sessions with su but no SC
% TODO: find long chains? loops?

critset=categorical({'Naive','Learning','WT'});
prefix=["nv","ln","wt"];
SAMP_SU=15;
MIN_SESS_CNT=20;
PER_SESS_RPT=100;

nlw_data.wt_sel_meta=wt_sel_meta;
nlw_data.nv_sel_meta=nv_sel_meta;
nlw_data.ln_sel_meta=ln_sel_meta;

nlw_data.wt_su_meta=wt_su_meta;
nlw_data.nv_su_meta=nv_su_meta;
nlw_data.ln_su_meta=ln_su_meta;

nlw_data.wt_sig=wt_sig;
nlw_data.nv_sig=nv_sig;
nlw_data.ln_sig=ln_sig;

outtbl=cell(0);

for crit_idx=1:3
    usess=intersect(unique(nlw_data.(prefix(crit_idx)+"_su_meta").sess), unique(nlw_data.(prefix(crit_idx)+"_sig").sess));
    for sidx=reshape(usess,1,[])
        if rem(sidx,10)==0
            disp([critset(crit_idx),num2str(sidx)]);
        end
        sess_sel=nlw_data.(prefix(crit_idx)+"_su_meta").sess==sidx;

        s1sel=sess_sel & nlw_data.(prefix(crit_idx)+"_sel_meta").wave_id==5;
        s2sel=sess_sel & nlw_data.(prefix(crit_idx)+"_sel_meta").wave_id==6;
        nmsel=sess_sel & nlw_data.(prefix(crit_idx)+"_sel_meta").wave_id==0;

        if all([nnz(s1sel),nnz(s2sel),nnz(nmsel)]>MIN_SESS_CNT)
            for rpt=1:PER_SESS_RPT
                s1_su=randsample(nlw_data.(prefix(crit_idx)+"_su_meta").allcid(s1sel),SAMP_SU,false);
                s2_su=randsample(nlw_data.(prefix(crit_idx)+"_su_meta").allcid(s2sel),SAMP_SU,false);
                nm_su=randsample(nlw_data.(prefix(crit_idx)+"_su_meta").allcid(nmsel),SAMP_SU,false);

                mem_su=[s1_su;s2_su];
                all_su=[mem_su;nm_su];

                sig_sel=nlw_data.(prefix(crit_idx)+"_sig").sess==sidx;
                sesssig=nlw_data.(prefix(crit_idx)+"_sig").suid(sig_sel,:);

                mem_su_sig_sel=all(ismember(sesssig,mem_su),2);
                mem_gh=digraph(string(sesssig(mem_su_sig_sel,1)),string(sesssig(mem_su_sig_sel,2)));

                s1_su_sig_sel=all(ismember(sesssig,s1_su),2);
                s1_gh=digraph(string(sesssig(s1_su_sig_sel,1)),string(sesssig(s1_su_sig_sel,2)));

                s2_su_sig_sel=all(ismember(sesssig,s2_su),2);
                s2_gh=digraph(string(sesssig(s2_su_sig_sel,1)),string(sesssig(s2_su_sig_sel,2)));

                nm_su_sig_sel=all(ismember(sesssig,nm_su),2);
                nm_gh=digraph(string(sesssig(nm_su_sig_sel,1)),string(sesssig(nm_su_sig_sel,2)));


                [~,s1_binsize]=s1_gh.conncomp('Type','weak');
                [~,s2_binsize]=s2_gh.conncomp('Type','weak');
                [~,mem_binsize]=mem_gh.conncomp('Type','weak');
                [~,nm_binsize]=nm_gh.conncomp('Type','weak');

                % starts with degrees and chain lengths
                currcell={critset(crit_idx),sidx,rpt,...
                    s1_binsize,s1_gh.indegree,s1_gh.outdegree,... % 4-6
                    s2_binsize,s2_gh.indegree,s2_gh.outdegree,... % 7-9
                    mem_binsize,mem_gh.indegree,mem_gh.outdegree,... % 10-12
                    nm_binsize,nm_gh.indegree,nm_gh.outdegree,...% 13-15
                    }; 

                outtbl=[outtbl;currcell];
            end
        end
    end
end


% per-component-size
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=[outtbl{ii,4},outtbl{ii,7}];
    kwv=[kwv;ss(ss>2).'];
    kwg=[kwg;repmat(outtbl{ii,1},nnz(ss>2),1)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
figure()
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Connected component size');
title("Congruent ensembles","p = "+num2str(p,4))


% per-component-size
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=outtbl{ii,10};
    kwv=[kwv;ss(ss>2).'];
    kwg=[kwg;repmat(outtbl{ii,1},nnz(ss>2),1)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
figure()
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Connected component size');
title("Ensembles of all memory neuron","p = "+num2str(p,4))

% per-component-sizekwv=[];
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=outtbl{ii,13};
    kwv=[kwv;ss(ss>2).'];
    kwg=[kwg;repmat(outtbl{ii,1},nnz(ss>2),1)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
figure()
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Connected component size');
title("Non-memory neuron","p = "+num2str(p,4))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out-degree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% per-component-size
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=[outtbl{ii,5};outtbl{ii,8}];
    kwv=[kwv;ss(ss>0)];
    kwg=[kwg;repmat(outtbl{ii,1},nnz(ss>0),1)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
figure()
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('In degree');
title("Congruent ensembles","p = "+num2str(p,4))


% per-component-size
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=outtbl{ii,10};
    kwv=[kwv;ss(ss>2).'];
    kwg=[kwg;repmat(outtbl{ii,1},nnz(ss>2),1)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
figure()
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Connected component size');
title("Ensembles of all memory neuron","p = "+num2str(p,4))

% per-component-sizekwv=[];
kwv=[];
kwg=[];
for ii=1:height(outtbl)
    ss=outtbl{ii,13};
    kwv=[kwv;ss(ss>2).'];
    kwg=[kwg;repmat(outtbl{ii,1},nnz(ss>2),1)];
end
% kruskalwallis(kwv,kwg)
[gidx,gid]=findgroups(kwg);
mm=splitapply(@(x) mean(x),kwv,gidx);
sem=splitapply(@(x) std(x)/sqrt(numel(x)),kwv,gidx);
figure()
bh=bar(diag(flip(mm)),'stacked');
hold on
eh=errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.');
set(gca(),'XTickLabel',{'Naive','Learning','Well-trained'})
p=anovan(kwv,{kwg},'display','off');
ylabel('Connected component size');
title("Non-memory neuron","p = "+num2str(p,4))




