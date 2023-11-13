flall=dir('K:\LiukaiData\**\*_offset.mat');
dl=unique({flall.folder}).';
session_su=[];
for di=1:numel(dl)
    fl=dir(fullfile(dl{di},'*_offset.mat'));
    spkid=[];
    for fi=1:numel(fl)
        tfstr=load(fullfile(fl(fi).folder,fl(fi).name));
        susel=arrayfun(@(x) ~isempty(tfstr.firings_offset(x).unit_indx),1:numel(tfstr.firings_offset));
        if ~any(susel)
            continue
        end
        for leadii=reshape(find(susel),1,[])
            spkid=[spkid;reshape(tfstr.firings_offset(leadii).unit_indx+1000*leadii+1000000*fi,[],1)];
        end
    end
    session_su=[session_su;cell2table({di,dl{di},spkid},'VariableNames',{'idx','path','suids'})];
end
total_su=cellfun(@(x) numel(x),session_su.suids);
session_su.total_su=total_su;



out=[];
%% load data
for di=89%1:numel(dl)
    fl=dir(fullfile(dl{di},'*_offset.mat'));
    spkts=[];
    spkid=[];
    for fi=1:numel(fl)
        tfstr=load(fullfile(fl(fi).folder,fl(fi).name));
        susel=arrayfun(@(x) ~isempty(tfstr.firings_offset(x).unit_indx),1:numel(tfstr.firings_offset));
        if ~any(susel)
            continue
        end
        for leadii=reshape(find(susel),1,[])
            nsu=numel(tfstr.firings_offset(leadii).unit_indx);
            for uii=1:nsu
                ts=tfstr.firings_offset(leadii).ts{uii}.';
                spkts=[spkts;double(ts)];
                id=tfstr.firings_offset(leadii).unit_indx(uii)+1000*leadii+1000000*fi;
                spkid=[spkid;double(repmat(id,numel(ts),1))];
                % spkid=[spkid;id];
            end
        end
    end

    % out=[out;cell2table({di,dl{di},spkid},'VariableNames',{'idx','path','suids'})];

end
%% FC
ephys.util.dependency()
rez=bz.sortSpikeIDz(spkts,spkid);

%% FR
beh_mat=readmatrix("K:\LiukaiData\20230908\01-wam-8-moving\01-wam-8-moving-ag-20230908_104404-2\01-wam-8-moving-ag-20230908_104404-2.DIO_din.csv");
stimOnsetSel=ismember(beh_mat(:,2),21:28) & beh_mat(:,3)==1;
onsetTS=beh_mat(stimOnsetSel,[1,1,2]);

% disp(mink(diff(onsetTS),10));% 20tick / ms
onsetTS(:,2)=onsetTS(:,1)+900*20;
edges=reshape(onsetTS(:,1:2).',[],1);
[~,~,trledge]=histcounts(spkts,edges);

%[trl#,loc,spkid,fr]
frmat=[];
for oneid=reshape(unique(spkid),1,[])
    for onetrl=1:size(onsetTS,1) % TODO: vectorize for performance
        oneloc=onsetTS(onetrl,3);
        fr=nnz(spkid==oneid & (trledge==onetrl*2-1))./0.9;
        frmat=[frmat;onetrl,oneloc,oneid,fr];
    end
end

%% selectivity
selmat=[];
for oneid=reshape(unique(frmat(:,3)),1,[])
    susel=frmat(:,3)==oneid;
    p=anova1(frmat(susel,4),frmat(susel,2),'off');
    selmat=[selmat;oneid,p];
end

