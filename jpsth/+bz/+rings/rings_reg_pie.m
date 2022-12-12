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

% olf pie
olfsel=strcmp(rstats(:,9),'olf');
tot=cell(0);
for ii=reshape(find(olfsel),1,[])
    tot=[tot,rstats{ii,8}];
end

unique(tot)




