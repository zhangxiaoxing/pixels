% find bridge neurons, aka chain neurons between loops
global_init;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);

%% single spike chain
sschain=load('chain_tag.mat','out');
keys=[struct2cell(structfun(@(x) fieldnames(x), sschain.out.d6, 'UniformOutput', false));...
    struct2cell(structfun(@(x) fieldnames(x), sschain.out.d3, 'UniformOutput', false))];
keys=vertcat(keys{:});
ssc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));

%% single spike loop
load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats'); % bz.rings.rings_time_constant
pstats=rmfield(pstats,"nonmem");
loopfn=fieldnames(pstats.congru);
loopsess=str2double(regexp(loopfn,'(?<=s)\d{1,3}(?=r)','match','once'));
ssl_sess=unique(loopsess);

usess=intersect(ssc_sess,ssl_sess);

% single spk chn:1, burst spk chn:2, single spk loop:4, burst spk loop:8

bridge_sums=[];
pivot_sums=[];
for sessid=reshape(usess,1,[])
   
    % all loop cid
    %% single spike loop
    loopset=[];
    for cc=reshape(loopfn(loopsess==sessid),1,[])
        loopset=[loopset,pstats.congru.(cc{1}).rstats{1,3}];
    end
    loopset=unique(loopset);

    %% single spike chain
    for dur=["d6","d3"]
        for wid=reshape(fieldnames(sschain.out.(dur)),1,[])
            for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[])
                if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                    continue
                end
                onechain=sschain.out.(dur).(wid{1}).(cc{1}).meta{1};
                for suidx=2:(numel(onechain)-1)
                    if any(ismember(onechain(1:suidx-1),loopset),"all") ...
                            && any(ismember(onechain(suidx+1:end),loopset),"all")...
                            && ~ismember(onechain(suidx),loopset)
                        bridge_sums=[bridge_sums;sessid,onechain(suidx)];
                        % pivot
                        if ismember(onechain(suidx-1),loopset)
                            pivot_sums=[pivot_sums;sessid,onechain(suidx-1)];
                        end
                        if ismember(onechain(suidx+1),loopset)
                            pivot_sums=[pivot_sums;sessid,onechain(suidx+1)];
                        end

                    end
                end
            end
        end
    end
end
bucid=unique(bridge_sums,"rows");
bregs=cell(size(bucid,1),1);
for ii=1:size(bucid,1)
    susel=su_meta.sess==bucid(ii,1) & su_meta.allcid==bucid(ii,2);
    bregs(ii)=su_meta.reg_tree(5,susel);
end
brregs=categorical(bregs(~ismissing(bregs)));

cucid=unique(pivot_sums,"rows");
cregs=cell(size(cucid,1),1);
for ii=1:size(cucid,1)
    susel=su_meta.sess==cucid(ii,1) & su_meta.allcid==cucid(ii,2);
    cregs(ii)=su_meta.reg_tree(5,susel);
end
crregs=categorical(cregs(~ismissing(cregs)));

figure()
tiledlayout(1,2);
nexttile()
pie(brregs)
title('Bridge type')
nexttile()
pie(crregs)
title('Anchor type')


