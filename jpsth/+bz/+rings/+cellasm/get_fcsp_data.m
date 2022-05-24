function out_=get_fcsp_data()
persistent out
if isempty(out)
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
    ctxset=[];
    for keyc=idmap.reg2tree.keys
        key=keyc{1};
        tree=idmap.reg2tree(key);
        if numel(tree)>6 && strcmp(tree{4},'CTX') && ~isempty(tree{7})
            ctxset=[ctxset;idmap.reg2ccfid(tree{7})];
        end
    end
    ctxset=unique(ctxset);
    [metas,~,~,per_trial]=bz.fccoding.get_fc_coding('per_trial',true);

    %meta:sess,suid1,suid2,mem,reg
    ctx_sel=all(ismember(metas(:,6:7),ctxset),2);
    s1congru=all(ismember(metas(:,4:5),1:2),2);
    s2congru=all(ismember(metas(:,4:5),3:4),2);

    out=cell(0,1);

    s1idces=find(ctx_sel & s1congru).';
    for fcid=s1idces
        unit_map=containers.Map('KeyType','char','ValueType','any');
        unit_map('pref')=per_trial{fcid,1};
        unit_map('nonp')=per_trial{fcid,2};
        unit_map('eprf')=per_trial{fcid,3};
        unit_map('enon')=per_trial{fcid,4};
        unit_map('meta')={'FCSP',2,fcid,metas(fcid,1),metas(fcid,2:3)};
        out{end+1,1}=unit_map;
    end

    s2idces=find(ctx_sel & s2congru).';
    for fcid=s2idces
        unit_map=containers.Map('KeyType','char','ValueType','any');
        unit_map('pref')=per_trial{fcid,2};
        unit_map('nonp')=per_trial{fcid,1};
        unit_map('eprf')=per_trial{fcid,4};
        unit_map('enon')=per_trial{fcid,3};
        unit_map('meta')={'FCSP',2,fcid,metas(fcid,1),metas(fcid,2:3)};
        out{end+1,1}=unit_map;
    end
end
out_=out;
end


