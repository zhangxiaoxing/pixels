% tag spikes in neuron with loop or chain activity
denovo=true;
skipfile=true;

if denovo
    inited=true;
    %% single spike chain
    sschain=load('chain_tag.mat','out');
    keys=[struct2cell(structfun(@(x) fieldnames(x), sschain.out.d6, 'UniformOutput', false));...
        struct2cell(structfun(@(x) fieldnames(x), sschain.out.d3, 'UniformOutput', false))];
    keys=vertcat(keys{:});
    ssc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));
    
    %% single spike loop
    load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats'); % bz.rings.rings_time_constant
    pstats=rmfield(pstats,"nonmem");
    ssl_sess=unique(str2double(regexp(fieldnames(pstats.congru),'(?<=s)\d{1,3}(?=r)','match','once')));
    usess=intersect(ssc_sess,ssl_sess);
    
    %% per-session entrance
    cid_motif_maps=struct();

    for sessid=reshape(usess,1,[])
        cid_motif_maps.("S"+sessid)=containers.Map('KeyType','int32','ValueType','any');
        [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
        wtsel=trials(:,9)>0 & trials(:,10)>0 & ismember(trials(:,5),[4 8]) &ismember(trials(:,8),[3 6]);

        %% single spike chain
        for dur=["d6","d3"]
            dd=str2double(replace(dur,"d",""));
            for wid=reshape(fieldnames(sschain.out.(dur)),1,[])
                for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[])
                    if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                        continue
                    end
                    onechain=sschain.out.(dur).(wid{1}).(cc{1});
                    for mcid=onechain.meta{1}
                        if ~cid_motif_maps.("S"+sessid).isKey(int32(mcid))
                            cid_motif_maps.("S"+sessid)(int32(mcid))=cell(0);
                        end
                        cid_motif_maps.("S"+sessid)(int32(mcid))=[cid_motif_maps.("S"+sessid)(int32(mcid));cc{1}];
                    end
                end
            end
        end
        %% single spike loop
        for cc=reshape(fieldnames(pstats.congru),1,[])
            if ~startsWith(cc{1},['s',num2str(sessid),'r'])
                continue
            end
            oneloop=pstats.congru.(cc{1});
            for mcid=oneloop.rstats{3}
                if ~cid_motif_maps.("S"+sessid).isKey(int32(mcid))
                    cid_motif_maps.("S"+sessid)(int32(mcid))=cell(0);
                end
                cid_motif_maps.("S"+sessid)(int32(mcid))=[cid_motif_maps.("S"+sessid)(int32(mcid));cc{1}];
            end
        end
    end
    if ~skipfile
        blame=vcs.blame();
        save(fullfile("bzdata","single_su_multi_motif.mat"),'cid_motif_maps','blame')
    end
    denovoflag=true;
else % load from file
    load(fullfile("bzdata","single_su_multi_motif.mat"),'cid_motif_maps','blame')
end

sums=[];
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
for ss=reshape(struct2cell(cid_motif_maps),1,[])
    sums=[sums,cellfun(@(x) numel(unique(x)), ss{1}.values)];
end
hhist=histcounts(sums,[0.5:600.5],'Normalization','cdf');
figure()
hold on;
plot(1:600,hhist,'-k');
% set(gca(),'YScale','log','XTick',0:50:150)
ylim([0,1])
xlim([0,600])
xlabel('Number of composite loops')
ylabel('Probability density')
title('Single SU in composite loops')
set(gca,'XScale','log')


%% previous implementation