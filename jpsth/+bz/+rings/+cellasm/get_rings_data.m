function out=get_rings_data()

persistent loops_sums meta
if isempty(loops_sums) || isempty(meta)
    load('loops_coding.mat','loops_sums')  %from \+bz\+rings\rings_coding.m
    meta=ephys.util.load_meta();
end
out=cell(0,1);
mtypes=["congru"];%["congru","nonmem"];
regtypes=["cross","within"];
rsizes=3:5;
for rsize=rsizes
    for mtype=mtypes
        for rtype=regtypes
            rings=cell2mat(loops_sums.(mtype).(sprintf('%s_%d',rtype,rsize)).');
            for ridx=1:numel(rings)
                if rem(ridx,1000)==0, disp([rsize,ridx]);end
                unit_map=containers.Map('KeyType','char','ValueType','any');
                trials=rings(ridx).trialinfo;
                pref_trials=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,9)==1 & trials(:,10)==1);
                nonp_trials=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,9)==1 & trials(:,10)==1);
                pref_err_trls=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,10)==0);
                nonp_err_trls=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,10)==0);
                %TODO room for improvment
                memtype=meta.mem_type(meta.sess==rings(ridx).sessid & meta.allcid==str2double(rings(ridx).label{1}));
                
                if strcmp(mtype,'congru') && memtype>2
                    [pref_trials,nonp_trials]=deal(nonp_trials,pref_trials); %S1 S2 swap
                    [pref_err_trls,nonp_err_trls]=deal(nonp_err_trls,pref_err_trls);
                end

                for ci=1:numel(rings(ridx).label)
                    unit_map('pref')=arrayfun(@(x) nnz(rings(ridx).trial{ci}==x & rings(ridx).time{ci} >= 1 & rings(ridx).time{ci} < 7),pref_trials);
                    unit_map('nonp')=arrayfun(@(x) nnz(rings(ridx).trial{ci}==x & rings(ridx).time{ci} >= 1 & rings(ridx).time{ci} < 7),nonp_trials);
                    unit_map('eprf')=arrayfun(@(x) nnz(rings(ridx).trial{ci}==x & rings(ridx).time{ci} >= 1 & rings(ridx).time{ci} < 7),pref_err_trls);
                    unit_map('enon')=arrayfun(@(x) nnz(rings(ridx).trial{ci}==x & rings(ridx).time{ci} >= 1 & rings(ridx).time{ci} < 7),nonp_err_trls);
                    unit_map('meta')={rtype,rsize,ridx,rings(ridx).sessid,str2double(rings(ridx).label)};
                    out{end+1,1}=unit_map;
                end
            end
        end
    end
end
end
