function loop_occurance_per_reg_su(sums_all,su_meta,opt)
arguments
    sums_all
    su_meta
    opt.pie (1,1) logical = true
    opt.bar (1,1) logical = true
    opt.loopPerSU (1,1) logical = true
end

% su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
% load(fullfile('bzdata','sums_ring_stats_all.mat'));

% greyregs=categorical(ephys.getGreyRegs('range','grey'));
su_reg=categorical(su_meta.reg_tree(5,:));
rstats=bz.rings.rings_reg_pie(sums_all,'plot',false);

olfsel=strcmp(rstats(:,9),'olf');
bothsel=strcmp(rstats(:,9),'both');
within_sel=cellfun(@(x) numel(unique(x)), rstats(:,8))==1;

loop_reg=categorical(unique([rstats{:,8}]));

olfreg_within=categorical([rstats{olfsel & within_sel,8}]);
olfreg_cross=categorical([rstats{olfsel & ~within_sel,8}]);
bothreg_within=categorical([rstats{bothsel & within_sel,8}]);
bothreg_cross=categorical([rstats{bothsel & ~within_sel,8}]);

id_olf_within=[rstats{olfsel & within_sel,10}];
[uid_olf_within,ia,~]=unique(id_olf_within);
uid_olf_within_reg=olfreg_within(ia);

id_olf_cross=[rstats{olfsel & ~within_sel,10}];
[uid_olf_cross,ia,~]=unique(id_olf_cross);
uid_olf_cross_reg=olfreg_cross(ia);

id_both_within=[rstats{bothsel & within_sel,10}];
[uid_both_within,ia,~]=unique(id_both_within);
uid_both_within_reg=bothreg_within(ia);

id_both_cross=[rstats{bothsel & ~within_sel,10}];
[uid_both_cross,ia,~]=unique(id_both_cross);
uid_both_cross_reg=bothreg_cross(ia);


lbls={{'Olfactory within region','Sum of duplicates in loops'};...
    {'Olfactory cross region','Sum of duplicates in loops'};...
    {'Olfactory within region','Unique neurons in loops'};...
    {'Olfactory cross region','Unique neurons in loops'};...
    {'Both within region','Sum of duplicates in loops'};...
    {'Both cross region','Sum of duplicates in loops'};...
    {'Both within region','Unique neurons in loops'};...
    {'Both cross region','Unique neurons in loops'}};
dsets={olfreg_within,olfreg_cross,uid_olf_within_reg,uid_olf_cross_reg,...
        bothreg_within,bothreg_cross,uid_both_within_reg,uid_both_cross_reg};

if opt.pie
    figure()
    tiledlayout(2,4)
    cnt=1;
    for dset=dsets
        nexttile()
        pie(dset{1})
        title(lbls{cnt});
        cnt=cnt+1;
    end
end

if opt.bar
    figure()
    tiledlayout(2,4)
    cnt=1;

    for dset=dsets
        nexttile()
        ratios=[];
        for creg=loop_reg
            reg_su_cnt=nnz(su_reg==creg);
            dset_cnt=nnz(dset{1}==creg);
            ratios=[ratios,dset_cnt./reg_su_cnt];
        end
        [sratios,sidx]=sort(ratios,'descend','MissingPlacement','last');
        if any(contains(lbls{cnt},'within'),'all')
            bar(sratios(1:2))
            set(gca(),"XTick",1:2,"XTickLabel",loop_reg(sidx(1:2)),...
            "XTickLabelRotation",90,'YScale','log','YLim',[0.001,1])
        else
            bar(sratios(1:5))
            set(gca(),"XTick",1:5,"XTickLabel",loop_reg(sidx(1:5)),...
            "XTickLabelRotation",90,'YScale','log','YLim',[0.001,1])
        end
        xlim([0.25,5.75]);

        ylabel('Occurance in loops per neuron')
        title(lbls{cnt});


        cnt=cnt+1;
    end
end

if opt.loopPerSU  % loops per su, classed by region
ddsets={id_olf_within,uid_olf_within,uid_olf_within_reg;...
    id_olf_cross,uid_olf_cross,uid_olf_cross_reg;...
       id_both_within,uid_both_within,uid_both_within_reg;...
       id_both_cross,uid_both_cross,uid_both_cross_reg};
lbls={{'Olfactory within region'};...
    {'Olfactory cross region'};...
    {'Both within region'};...
    {'Both cross region'}};

    
    
    figure()
    tiledlayout(2,2)
    for didx=1:4
        nexttile();
        mean_occurance=[];
        for creg=loop_reg
            reg_sel=ddsets{didx,3}==creg;
            if nnz(reg_sel)>0
                uids=ddsets{didx,2}(reg_sel);
                loops_cnt=nnz(ismember(ddsets{didx,1},uids));
                mean_occurance=[mean_occurance;loops_cnt./nnz(reg_sel)];
            else
                mean_occurance=[mean_occurance;0];
            end
        end
        bar(mean_occurance)
        set(gca(),"XTick",1:11,"XTickLabel",loop_reg,"XTickLabelRotation",90)
        ylabel('Occurance in loops per unique neuron')
        title(lbls{didx});
    end
end

end