function out=get_SU_data()

persistent meta
if isempty(meta)
    meta=ephys.util.load_meta();
end
sens_sel=meta.mem_type>0;
ctx_sel=strcmp(meta.reg_tree(2,:),'CTX') & cellfun(@(x) length(x),meta.reg_tree(5,:))>0;
global_sel=(sens_sel & ctx_sel).';
sesses=unique(meta.sess(global_sel));

out=cell(0,1);
homedir=ephys.util.getHomedir('type','raw');
[~,~,sessmap]=ephys.sessid2path(0);
for ii=reshape(sesses,1,[])
    disp(ii)
    fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
    fr=h5read(fpath,'/FR_All');
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    susel=global_sel(meta.sess==ii);
    memtypes=meta.mem_type(meta.sess==ii);

    fr=fr(:,susel,:);
    suid=suid(susel);
    memtypes=memtypes(susel);

    pref_trials=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,9)==1 & trials(:,10)==1);
    nonp_trials=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,9)==1 & trials(:,10)==1);
    pref_err_trls=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,10)==0);
    nonp_err_trls=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,10)==0);
    
    for ci=1:size(fr,2)
        if memtypes(ci)>2
            [pref_trials,nonp_trials]=deal(nonp_trials,pref_trials); %S1 S2 swap
            [pref_err_trls,nonp_err_trls]=deal(nonp_err_trls,pref_err_trls);
        end
        unit_map=containers.Map('KeyType','char','ValueType','any');
        unit_map('pref')=arrayfun(@(x) mean(fr(x,ci,17:40),3),pref_trials);
        unit_map('nonp')=arrayfun(@(x) mean(fr(x,ci,17:40),3),nonp_trials);
        unit_map('eprf')=arrayfun(@(x) mean(fr(x,ci,17:40),3),pref_err_trls);
        unit_map('enon')=arrayfun(@(x) mean(fr(x,ci,17:40),3),nonp_err_trls);
        unit_map('meta')={'su',1,ci,ii,suid(ci)};
        out{end+1,1}=unit_map;
    end
end

end
