su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
load(fullfile('bzdata','sums_ring_stats_all.mat'));
rstats=cell(0,10);
for rsize=3:5
    one_rsize=sums_all{rsize-2};
    curr_sess=-1;
    for ridx=1:size(one_rsize,1)
        if curr_sess~=one_rsize{ridx,1}
            curr_sess=one_rsize{ridx,1};
            sesscid=su_meta.allcid(su_meta.sess==curr_sess);
            sesswaveid=wrs_mux_meta.wave_id(su_meta.sess==curr_sess);
            sess_part=su_meta.reg_tree(1,su_meta.sess==curr_sess);
            sess_reg=su_meta.reg_tree(5,su_meta.sess==curr_sess);

            sess_wave_map=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
            sess_part_map=containers.Map(num2cell(sesscid),sess_part);
            sess_reg_map=containers.Map(num2cell(sesscid),sess_reg);
        end
        curr_waveid=cell2mat(sess_wave_map.values(num2cell(one_rsize{ridx,3})));
        curr_part=sess_part_map.values(num2cell(one_rsize{ridx,3}));
        curr_reg=sess_reg_map.values(num2cell(one_rsize{ridx,3}));
        if (~all(ismember(curr_part,{'CH','BS'}),'all')) || any(ismissing(curr_reg),'all')
            continue
        end
        [rwid,seltype]=bz.rings.ring_wave_type(curr_waveid);
        if strcmp(rwid,'congru')
            rstats=[rstats;one_rsize(ridx,:),curr_waveid,{curr_reg},seltype];
        end
    end
end

rstats(:,10)=arrayfun(@(x) rstats{x,1}*100000+rstats{x,3},1:size(rstats,1),'UniformOutput',false);

% olf pie
olfsel=strcmp(rstats(:,9),'olf');
bothsel=strcmp(rstats(:,9),'both');
olfreg=categorical([rstats{olfsel,8}]);
bothreg=categorical([rstats{bothsel,8}]);

uid=[rstats{olfsel,10}];
[~,ia,~]=unique(uid);
uid_olf_reg=olfreg(ia);

uid=[rstats{bothsel,10}];
[~,ia,~]=unique(uid);
uid_both_reg=bothreg(ia);



figure()
tiledlayout(2,2)
nexttile()
ph=pie(olfreg);
for hh=reshape(ph,1,[])
    if isprop(hh,'EdgeColor')
        hh.EdgeColor='none';
    end
end
legend(categories(olfreg));
title({'Olfactory','Sum of duplicates in loops'})

nexttile()
ph=pie(uid_olf_reg);
for hh=reshape(ph,1,[])
    if isprop(hh,'EdgeColor')
        hh.EdgeColor='none';
    end
end
title({'Olfactory','Unique neurons in loops'})

nexttile()
ph=pie(bothreg);
for hh=reshape(ph,1,[])
    if isprop(hh,'EdgeColor')
        hh.EdgeColor='none';
    end
end
legend(categories(bothreg));
title({'Both-selective','Sum of duplicates in loops'})

nexttile()
ph=pie(uid_both_reg);
for hh=reshape(ph,1,[])
    if isprop(hh,'EdgeColor')
        hh.EdgeColor='none';
    end
end
title({'Both-selective','Unique neurons in loops'})

function total_SU_count

end
