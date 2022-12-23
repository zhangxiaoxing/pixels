% chains_uf=wave.COM_chain(sel_meta);
% chains_uf_rev=wave.COM_chain(sel_meta,'reverse',true);
% blame=vcs.blame();
% save('chains_mix.mat','chains_uf','chains_uf_rev','blame')

global_init;
load('chains_mix.mat','chains_uf','chains_uf_rev');

% vs shuffle
% jpsth/+bz/+rings/shuffle_conn_bz_alt.m
% wave.COM_chain_shuf(wrs_mux_meta);
load('chains_shuf.mat','shuf_chains')

su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
chains_uf.reg=cell(size(chains_uf.cids));
chains_uf_rev.reg=cell(size(chains_uf_rev.cids));
chains_uf.reg_sel=false(size(chains_uf.cids));
chains_uf_rev.reg_sel=false(size(chains_uf_rev.cids));
greys=ephys.getGreyRegs('range','grey');

for sess=reshape(unique(chains_uf.sess),1,[])
    sesscid=su_meta.allcid(su_meta.sess==sess);
    sessreg=su_meta.reg_tree(5,su_meta.sess==sess);
    for cnid=reshape(find(chains_uf.sess==sess),1,[])
        [~,supos]=ismember(chains_uf.cids{cnid},sesscid);
        chains_uf.reg{cnid}=sessreg(supos);
        if all(ismember(sessreg(supos),greys),'all')
            chains_uf.reg_sel(cnid)=true;
        end
    end
    for cnid=reshape(find(chains_uf_rev.sess==sess),1,[])
        [~,supos]=ismember(chains_uf_rev.cids{cnid},sesscid);
        chains_uf_rev.reg{cnid}=sessreg(supos);
        if all(ismember(sessreg(supos),greys),'all')
            chains_uf_rev.reg_sel(cnid)=true;
        end
    end
end %sess

%% total number
countbin=2.5:12.5;
data_hist=histcounts(cellfun(@(x) numel(x),chains_uf.cids(chains_uf.reg_sel)),countbin);
data_rev_hist=histcounts(cellfun(@(x) numel(x),chains_uf_rev.cids(chains_uf_rev.reg_sel)),countbin);
if false
    load('chains_nonmem.mat','nonmem_chains')
    nonmem_hist=histcounts(cellfun(@(x) numel(x),nonmem_chains.cids),countbin);
end

shuf_hist=nan(numel(shuf_chains),numel(countbin)-1);
for shufid=1:numel(shuf_chains)
    shuf_hist(shufid,:)=histcounts(cellfun(@(x) numel(x),shuf_chains{shufid}.cids),countbin);    
end

shufmm=mean(shuf_hist);
shufsd=std(shuf_hist);

pxx=3:12;
figure()
hold on
fill([pxx,fliplr(pxx)],[shufmm-2*shufsd,fliplr(shufmm+2*shufsd)],'k','EdgeColor','none','FaceAlpha',0.2);
datah=plot(pxx,data_hist,'-r');
revh=plot(pxx,data_rev_hist,'-b');
shufh=plot(pxx,shufmm,'-k');
if true
    nonmemh=plot(pxx,nonmem_hist,'-','Color',[0.5,0.5,0.5]);
    legend([datah,revh,nonmemh,shufh],...
        {'Wave-consistent','Wave-inconsistent','Non-memory','Shuffled wave-consistent'});

else
    legend([datah,revh,shufh],{'Wave-consistent','Wave-inconsistent','Shuffled wave-consistent'});
end

xlabel('Number of linked neuron')
ylabel('Total occurance')
title('Total number of linked WM-neuron-series')

%% wave vs wave, 3s vs 6s
countbin=2.5:12.5;
%mix 3s, mix6s, olf 3s, olf 6s, dur 3s, dur 6s
wavetype={{'olf_s1','olf_s2'},{'dur_d3'},{'s1d3','s2d3'};...
    {'olf_s1','olf_s2'},{'dur_d6'},{'s1d6','s2d6'}};
durs=[3,3,3;6,6,6];
% ttls={'Olfactory, 3s','Duration, 3s','Both, 3s','Olfactory, 6s','Duration, 6s','Both, 6s'};
figure()
tiledlayout(1,3);
pxx=3:12;

lss={'--','-'};
for chartIdx=1:3
    nexttile;
    hold on
    for durIdx=1:2
        datasel=ismember(chains_uf.wave,wavetype{durIdx,chartIdx}) & chains_uf.dur==durs(durIdx,chartIdx);
        data_hist=histcounts(cellfun(@(x) numel(x),chains_uf.cids(chains_uf.reg_sel & datasel)),countbin);

        data_rev_sel=ismember(chains_uf_rev.wave,wavetype{durIdx,chartIdx}) & chains_uf_rev.dur==durs(durIdx,chartIdx);
        data_rev_hist=histcounts(cellfun(@(x) numel(x),chains_uf_rev.cids(chains_uf_rev.reg_sel & data_rev_sel)),countbin);

        shuf_hist=nan(numel(shuf_chains),numel(countbin)-1);
        for shufid=1:numel(shuf_chains)
            shuf_sel=ismember(shuf_chains{shufid}.wave,wavetype{durIdx,chartIdx}) ...
                & shuf_chains{shufid}.dur==durs(durIdx,chartIdx);
            shuf_hist(shufid,:)=histcounts(cellfun(@(x) numel(x),shuf_chains{shufid}.cids(shuf_sel)),countbin);
        end

        shufmm=mean(shuf_hist);
        shufsd=std(shuf_hist);

        fill([pxx,fliplr(pxx)],[shufmm-2*shufsd,fliplr(shufmm+2*shufsd)],'k','EdgeColor','none','FaceAlpha',0.2);

        datah=plot(pxx,data_hist,'r','LineStyle',lss{durIdx});
        revh=plot(pxx,data_rev_hist,'b','LineStyle',lss{durIdx});
        shufh=plot(pxx,shufmm,'k','LineStyle',lss{durIdx});
    end
    legend([datah,revh,shufh],{'Wave-consistent','Wave-inconsistent','Shuffled wave-consistent'});

    xlim([3,11])
    ylim([0,max(ylim())])
    xlabel('Number neuron in wave-link')
    ylabel('Total occurance')
%     title(ttls{chartIdx});
end


%% per-wave stats-per-region stats


% repeated occurance
chains_uf.uid=arrayfun(@(x) chains_uf.sess(x)*100000+int32(chains_uf.cids{x}), 1:numel(chains_uf.sess),'UniformOutput',false);

olf_sel=chains_uf.reg_sel & ismember(chains_uf.wave,["olf_s1","olf_s2"]);
olf_uid=[chains_uf.uid{olf_sel}];
[chain_olf_uid,usel]=unique(olf_uid);
olf_reg=[chains_uf.reg{olf_sel}];
olf_reg_cat=categorical(olf_reg);
olf_u_reg_cat=categorical(olf_reg(usel));


both_sel=chains_uf.reg_sel & ismember(chains_uf.wave,["s1d3","s2d3","s1d6","s2d6"]);
both_uid=[chains_uf.uid{both_sel}];
[chain_both_uid,usel]=unique(both_uid);
both_reg=[chains_uf.reg{both_sel}];
both_reg_cat=categorical(both_reg);
both_u_reg_cat=categorical(both_reg(usel));


dur_sel=chains_uf.reg_sel & ismember(chains_uf.wave,["dur_d3","dur_d6"]);
dur_uid=[chains_uf.uid{dur_sel}];
[chain_dur_uid,usel]=unique(dur_uid);
dur_reg=[chains_uf.reg{dur_sel}];
dur_reg_cat=categorical(dur_reg);
dur_u_reg_cat=categorical(dur_reg(usel));

figure()
tiledlayout(2,3)
nexttile()
olfh=pie(olf_reg_cat);
title("Olfactory with repeats")
nexttile()
durh=pie(dur_reg_cat);
title("Duration with repeats")
nexttile()
mixh=pie(both_reg_cat);
title("Both-selective with repeats")

nexttile()
olfh=pie(olf_u_reg_cat);
title("Olfactory unique neuron")
nexttile()
durh=pie(dur_u_reg_cat);
title("Duration unique neuron")
nexttile()
mixh=pie(both_u_reg_cat);
title("Both-selective unique neuron")


%chain vs loop
load(fullfile('bzdata','sums_ring_stats_all.mat'));
rstats=bz.rings.rings_reg_pie(sums_all,'plot',false);

lolf_sel=strcmp(rstats(:,9),'olf');
loop_olf_uid=unique([rstats{lolf_sel,10}]);
olf_shared=nnz(ismember(chain_olf_uid,loop_olf_uid));

lboth_sel=strcmp(rstats(:,9),'both');
loop_both_uid=unique([rstats{lboth_sel,10}]);
both_shared=nnz(ismember(chain_both_uid,loop_both_uid));

ldur_sel=strcmp(rstats(:,9),'dur');
loop_dur_uid=unique([rstats{ldur_sel,10}]);
dur_shared=nnz(ismember(chain_dur_uid,loop_dur_uid));

figure()
bh=bar([numel(chain_olf_uid)-olf_shared,olf_shared,numel(loop_olf_uid)-olf_shared;...
    numel(chain_dur_uid)-olf_shared,dur_shared,numel(loop_dur_uid)-dur_shared;...
    numel(chain_both_uid)-olf_shared,both_shared,numel(loop_both_uid)-both_shared],...
    'stacked');
bh(1).FaceColor='b';
bh(2).FaceColor='m';
bh(3).FaceColor='r';
legend(bh,{'Linked series only','Shared','Loops only'},'Location','northoutside','Orientation','horizontal','FontSize',12)
ylim([0,1800])
xlim([0.5,3.5])
ylabel('Number of neurons (out of 24667)')
set(gca(),'XTick',1:3,'XTickLabel',{'Olfactory','Duration','Both'},'XTickLabelRotation',30,'YTick',0:400:1600,'FontSize',12)

%% non-mem
% load('chains_nonmem.mat','nonmem_chains')
% countbin=2.5:12.5;
% data_hist=histcounts(cellfun(@(x) numel(x),nonmem_chains.cids),countbin);



%% TODO: performance-correlation
% TODO: per trial per spike align

%% TODO: wave time correlation?
chains_uf.reg_tcom=cell(size(chains_uf.reg));
sel6=chains_uf.dur==6 & chains_uf.reg_sel;
sel3=chains_uf.dur==3 & chains_uf.reg_sel;
olf_sel=ismember(chains_uf.wave,{'olf_s1','olf_s2'});
dur_sel=ismember(chains_uf.wave,{'dur_d3','dur_d6'});
both_sel=ismember(chains_uf.wave,{'s1d3','s2d3','s1d6','s2d6'});
%tcom_map:{mixed olf dur}

inkey=cellfun(@(x) all(ismember(x,tcom6_maps{2}.keys()),'all'),chains_uf.reg);
chains_uf.reg_tcom(sel6 & olf_sel & inkey)=cellfun(@(x) cell2mat(tcom6_maps{2}.values(x)),chains_uf.reg(sel6 & olf_sel & inkey),'UniformOutput',false);
xxs=cell2mat(cellfun(@(x) [x(1),x(end)],chains_uf.tcoms(sel6 & olf_sel & inkey),'UniformOutput',false));
yys=cell2mat(cellfun(@(x) [x(1),x(end)],chains_uf.reg_tcom(sel6 & olf_sel & inkey),'UniformOutput',false));
cross_sel=cellfun(@(x) numel(unique(x))>1,chains_uf.reg(sel6 & olf_sel & inkey));
long_sel=cellfun(@(x) numel(x)>4,chains_uf.reg(sel6 & olf_sel & inkey));

figure()
hold on
plot((xxs(cross_sel & long_sel,:).')./4,yys(cross_sel & long_sel,:).','-','Color',[0.49,0.5,0.51])
plot((xxs(cross_sel & long_sel,1).')./4,yys(cross_sel & long_sel,1).','o','MarkerFaceColor','b','MarkerEdgeColor','none')
plot((xxs(cross_sel & long_sel,2).')./4,yys(cross_sel & long_sel,2).','o','MarkerFaceColor','r','MarkerEdgeColor','none')

inkey=cellfun(@(x) all(ismember(x,tcom6_maps{3}.keys()),'all'),chains_uf.reg);
chains_uf.reg_tcom(sel6 & dur_sel & inkey)=cellfun(@(x) tcom6_maps{3}.values(x),chains_uf.reg(sel6 & dur_sel & inkey),'UniformOutput',false);

inkey=cellfun(@(x) all(ismember(x,tcom6_maps{1}.keys()),'all'),chains_uf.reg);
chains_uf.reg_tcom(sel6 & both_sel & inkey)=cellfun(@(x) tcom6_maps{1}.values(x),chains_uf.reg(sel6 & both_sel & inkey),'UniformOutput',false);


inkey=cellfun(@(x) all(ismember(x,tcom3_maps{2}.keys()),'all'),chains_uf.reg);
chains_uf.reg_tcom(sel3 & olf_sel & inkey)=cellfun(@(x) cell2mat(tcom3_maps{2}.values(x)),chains_uf.reg(sel3 & olf_sel & inkey),'UniformOutput',false);

xxs=cell2mat(cellfun(@(x) [x(1),x(end)],chains_uf.tcoms(sel3 & olf_sel & inkey),'UniformOutput',false));
yys=cell2mat(cellfun(@(x) [x(1),x(end)],chains_uf.reg_tcom(sel3 & olf_sel & inkey),'UniformOutput',false));
cross_sel=cellfun(@(x) numel(unique(x))>1,chains_uf.reg(sel3 & olf_sel & inkey));
long_sel=cellfun(@(x) numel(x)>4,chains_uf.reg(sel3 & olf_sel & inkey));

% figure()
% hold on
plot((xxs(cross_sel & long_sel,:).')./4,yys(cross_sel & long_sel,:).','--','Color',[0.49,0.5,0.51])
plot((xxs(cross_sel & long_sel,1).')./4,yys(cross_sel & long_sel,1).','o','MarkerFaceColor','b','MarkerEdgeColor','none')
plot((xxs(cross_sel & long_sel,2).')./4,yys(cross_sel & long_sel,2).','o','MarkerFaceColor','r','MarkerEdgeColor','none')



inkey=cellfun(@(x) all(ismember(x,tcom3_maps{3}.keys()),'all'),chains_uf.reg);
chains_uf.reg_tcom(sel3 & dur_sel & inkey)=cellfun(@(x) tcom3_maps{3}.values(x),chains_uf.reg(sel3 & dur_sel & inkey),'UniformOutput',false);

inkey=cellfun(@(x) all(ismember(x,tcom3_maps{1}.keys()),'all'),chains_uf.reg);
chains_uf.reg_tcom(sel3 & both_sel & inkey)=cellfun(@(x) tcom3_maps{1}.values(x),chains_uf.reg(sel3 & both_sel & inkey),'UniformOutput',false);






