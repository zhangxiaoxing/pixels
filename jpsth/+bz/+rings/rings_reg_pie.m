% TODO: pivot input to un-tagged format, for compatibility with shuffled
% data

function rstats=rings_reg_pie(rings,su_meta,opt)
arguments
    rings
    su_meta
    opt.plot (1,1) logical=true
    opt.congru (1,1) logical=true
end

wrs_mux_meta=ephys.get_wrs_mux_meta();
rstats=cell(0,10);
for rsize=3:5
    one_rsize=rings(:,rsize-2);
    for curr_sess=1:size(one_rsize,1)
        if isempty(rings{curr_sess,rsize-2})
            continue
        end
        sesscid=su_meta.allcid(su_meta.sess==curr_sess);
        sesswaveid=wrs_mux_meta.wave_id(su_meta.sess==curr_sess);
        sess_part=su_meta.reg_tree(1,su_meta.sess==curr_sess);
        sess_reg=su_meta.reg_tree(5,su_meta.sess==curr_sess);

        sess_wave_map=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
        sess_part_map=containers.Map(num2cell(sesscid),sess_part);
        sess_reg_map=containers.Map(num2cell(sesscid),sess_reg);
        for insessid=1:size(rings{curr_sess,rsize-2},1)
            curr_waveid=cell2mat(sess_wave_map.values(num2cell(rings{curr_sess,rsize-2}(insessid,:))));
            curr_part=sess_part_map.values(num2cell(rings{curr_sess,rsize-2}(insessid,:)));
            curr_reg=sess_reg_map.values(num2cell(rings{curr_sess,rsize-2}(insessid,:)));
            if (~all(ismember(curr_part,{'CH','BS'}),'all')) || any(ismissing(curr_reg),'all')
                continue
            end
            [rwid,seltype]=bz.rings.ring_wave_type(curr_waveid);
            if strcmp(rwid,'congru') || ~opt.congru
                rstats=[rstats;{rings{curr_sess,rsize-2}(insessid,:)},curr_waveid,{curr_reg},seltype,{curr_sess*100000+rings{curr_sess,rsize-2}(insessid,:)}];
            end
        end
    end
end


if opt.plot

    % olf pie
    olfsel=strcmp(rstats(:,4),'olf');
    bothsel=strcmp(rstats(:,4),'both');
    olfreg=categorical([rstats{olfsel,3}]);
    bothreg=categorical([rstats{bothsel,3}]);

    uid=[rstats{olfsel,5}];
    [~,ia,~]=unique(uid);
    uid_olf_reg=olfreg(ia);

    uid=[rstats{bothsel,5}];
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
end
end
% function total_SU_count
%
% end
