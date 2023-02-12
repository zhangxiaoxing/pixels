%% single spike chain
sschain=load('chain_tag.mat','out');
%% multi spike
bschain=load('chain_sust_tag_600.mat','out');
%% single spike loop
load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats');
pstats=rmfield(pstats,"nonmem");
%% burst spike loop, keys only
dbfile=fullfile("bzdata","rings_wave_burst_600.db");
% single spk chn:1, burst spk chn:2, single spk loop:4, burst spk loop:8

sessid=18;
load("ChainedLoop"+num2str(sessid)+".mat",'covered','FT_SPIKE')
edges = find(diff([0;covered;0]==1));
onset = edges(1:2:end-1);  % Start indices
run_length = edges(2:2:end)-onset;  % Consecutive ones counts

[~,maxid]=max(run_length);
maxonset=onset(maxid).*30;
maxoffset=edges(maxid*2).*30;
trialid=find(maxonset-FT_SPIKE.trialinfo(:,1)>0,1,'last');



[~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);

onsetTS=FT_SPIKE.trialinfo(trialid,1)+30000;
offsetTS=FT_SPIKE.trialinfo(trialid,1)+FT_SPIKE.trialinfo(trialid,8).*30000+30000;

% FT_SPIKE.lc_tag=cell(size(FT_SPIKE.timestamp));
% for cidx=1:numel(FT_SPIKE.timestamp)
%     FT_SPIKE.lc_tag{cidx}=zeros(size(FT_SPIKE.timestamp{cidx}));
% end
in_trial_tags=cell(0);
in_maximum_tags=cell(0);

%% single spike chain
for dur=["d6","d3"]
    for wid=reshape(fieldnames(sschain.out.(dur)),1,[])
        for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[])
            if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                continue
            end
            onechain=sschain.out.(dur).(wid{1}).(cc{1});
%             onechain=onechainstr.ts_id(onechainstr.ts_id(:,5)==trialid,:)
            occursel=find(all(onechain.ts>=onsetTS & onechain.ts<offsetTS,2));
            if ~isempty(occursel)
                maxsel=find(all(onechain.ts>=maxonset & onechain.ts<maxoffset));
                intrialsel=setdiff(occursel,maxsel);
%                 keyboard()
                if ~isempty(intrialsel)
                    for ii=reshape(intrialsel,1,[])
                        in_trial_tags(end+1,:)={'SSC',ii,[onechain.meta{1};onechain.ts(ii,:)].'};
                    end
                end
                if ~isempty(maxsel)
                    for ii=reshape(maxsel,1,[])
                        in_maximum_tags(end+1,:)={'SSC',ii,[onechain.meta{1};onechain.ts(ii,:)].'};
                    end
                end
            end
        end
    end
end

%% multi spike chain
for dur=["d6","d3"]
    for wid=reshape(fieldnames(bschain.out.(dur)),1,[])
        for cc=reshape(fieldnames(bschain.out.(dur).(wid{1})),1,[])
            if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                continue
            end
            onechain=bschain.out.(dur).(wid{1}).(cc{1});
            occursel=find(cellfun(@(x) all(x(:,3)>=onsetTS,'all') && all(x(:,3)<offsetTS,'all'),onechain.ts));
            if ~isempty(occursel)
                for ii=reshape(occursel,1,[])
                    if all(onechain.ts{ii}(:,3)>=maxonset,'all') && all(onechain.ts{ii}(:,3)<maxoffset,'all') 
                        in_maximum_tags(end+1,:)={'BSC',ii,...
                            [arrayfun(@(x) onechain.meta{1}(x), onechain.ts{ii}(:,1)),onechain.ts{ii}(:,3)]};
                    else
                        in_trial_tags(end+1,:)={'BSC',ii,...
                            [arrayfun(@(x) onechain.meta{1}(x), onechain.ts{ii}(:,1)),onechain.ts{ii}(:,3)]};
                    end
                end
            end
        end
    end
end


%% single spike loop
bracket=@(x) all(x>=onsetTS,'all') && all(x<offsetTS,'all');
for cc=reshape(fieldnames(pstats.congru),1,[])
    if ~startsWith(cc{1},['s',num2str(sessid),'r'])
        continue
    end
    onechain=pstats.congru.(cc{1});
    occursel=setdiff(arrayfun(@(x) x.*bracket(onechain.ts_id(onechain.ts_id(:,6)==x,1)),...
        setdiff(unique(onechain.ts_id(:,6)),0)),0);
    if ~isempty(occursel)
        for ii=reshape(occursel,1,[])
            curr=onechain.ts_id(onechain.ts_id(:,6)==ii,[2,1]);
            if all(curr(:,2)>=maxonset,'all') && all(curr(:,2)<maxoffset,'all')
                in_maximum_tags(end+1,:)={'SSL',ii,curr};
            else
                in_trial_tags(end+1,:)={'SSL',ii,curr};
            end
        end
    end
end

%% burst spike loop
dbfile=fullfile("bzdata","rings_wave_burst_600.db");
conn=sqlite(dbfile,"readonly");
bslkeys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
for cc=reshape(bslkeys,1,[])
    if ~contains(cc,['s',num2str(sessid),'r']) || ~endsWith(cc,'_ts')
        continue
    end
    chainmeta=table2array(conn.sqlread(replace(cc,'_ts','_meta')));
    maxrid=table2array(conn.fetch("SELECT MAX(Var1) from "+cc));
    for rid=1:1000:maxrid
        onechain=table2array(conn.fetch("SELECT * FROM "+cc+" WHERE Var1>="+num2str(rid)+ " AND Var1<"+num2str(rid+1000)));
        occursel=setdiff(arrayfun(@(x) x.*bracket(onechain(onechain(:,1)==x,4)),unique(onechain(:,1))),0);
        if ~isempty(occursel)
            for ii=reshape(occursel,1,[])
                curr=[arrayfun(@(x) chainmeta(x),onechain(onechain(:,1)==ii,2)),onechain(onechain(:,1)==ii,4)];
                if all(curr(:,2)>=maxonset,'all') && all(curr(:,2)<maxoffset,'all')
                    in_maximum_tags(end+1,:)={'BSL',ii,curr};
                else
                    in_trial_tags(end+1,:)={'BSL',ii,curr};
                end
            end
        end
    end
end

%% ===============PLOT PREPARATION===============

su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
sess_cid=su_meta.allcid(su_meta.sess==sessid);
sess_reg=su_meta.reg_tree(5,su_meta.sess==sessid).';

ucid=unique(subsref(cell2mat(in_maximum_tags(:,3)),struct(type={'()'},subs={{':',1}})));
[~,cidmap]=ismember(ucid,sess_cid);
ureg=sess_reg(cidmap);

ucid=ucid(~ismissing(ureg));
ureg=ureg(~ismissing(ureg));

% regional TCOM related 
global_init;
wrs_mux_meta=ephys.get_wrs_mux_meta();
if false
    com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false);
    [fcom6.olf.collection,fcom6.olf.com_meta]=wave.per_region_COM(...
        com_map,'sel_type','olf','com_field','com6');
    [~,tcidx]=ismember(ureg,fcom6.olf.collection(:,2));
    reg_tcom=cell2mat(fcom6.olf.collection(tcidx,1))./4;
    [reg_tcom,tcomidx]=sort(reg_tcom,'descend');
else
    [map_cells,~]=ephys.pct_reg_bars(wrs_mux_meta,'skip_plot',true); % only need map_cells for tcom-frac corr
    reg_prop=subsref(cell2mat(map_cells{2}.values(ureg)),struct(type={'()'},subs={{':',1}}));
    [reg_prop,tcomidx]=sort(reg_prop);
end

ureg=ureg(tcomidx);
ucid=ucid(tcomidx);

%% Actual plot

intrialmat=unique(cell2mat(in_trial_tags(:,3)),'rows');
inmaxmat=unique(cell2mat(in_maximum_tags(:,3)),'rows');

figure();
hold on;
yy=1;
cidyymap=containers.Map('KeyType','double','ValueType','double');
for cid=reshape(ucid,1,[])
    cidyymap(cid)=yy;
    ft_sel=strcmp(FT_SPIKE.label,num2str(cid));
    
    % TODO set diff tagged spike from all spikes
    ts=intrialmat(intrialmat(:,1)==cid,2);
    [~,tsidxx]=ismember(ts,FT_SPIKE.timestamp{ft_sel});
    trialxx=FT_SPIKE.time{ft_sel}(tsidxx);
    plot(trialxx,repmat(yy,1,numel(trialxx)),'|','MarkerEdgeColor','#FF80FF','LineWidth',1);

    ts=inmaxmat(inmaxmat(:,1)==cid,2);
    [~,tsidxx]=ismember(ts,FT_SPIKE.timestamp{ft_sel});
    maxxx=FT_SPIKE.time{ft_sel}(tsidxx);
    plot(maxxx,repmat(yy,1,numel(maxxx)),'|','MarkerEdgeColor','#FF0000','LineWidth',1);
    
    xx=FT_SPIKE.time{ft_sel}(FT_SPIKE.trial{ft_sel}==trialid);
    xx=setdiff(xx,union(trialxx,maxxx));
    plot(xx,repmat(yy,1,numel(xx)),'|','MarkerEdgeColor',"#808080")
    yy=yy+1;
end
% xlim([0,FT_SPIKE.trialinfo(trialid,8)]);
xlim([5,6.2]);
xlabel('Time (s)')
ylabel('Neuron #')
title("Session "+num2str(sessid)+" Trial "+num2str(trialid));
set(gca,'YTick',1:numel(ureg),'YTickLabel',ureg)

%% plot assembly
already=[];
cc=colormap('lines');
for ii=1:size(in_maximum_tags,1)
    curr=in_maximum_tags{ii,3};
    currfc=cell2mat(arrayfun(@(jj) [curr(jj,2),curr(jj,1),curr(jj+1,2),curr(jj+1,1)],(1:size(curr,1)-1).','UniformOutput',false));
    % skip unnecessary
    if ~isempty(already) && all(ismember(currfc,already,'rows'))
        disp("skipped "+num2str(ii))
        continue
    end
    already=unique([already;currfc],'rows');


    plotmat=[];
    for jj=1:size(curr,1)-1
        fromy=cidyymap(curr(jj,1));
        fromsel=strcmp(FT_SPIKE.label,num2str(curr(jj,1)));
        fromx=FT_SPIKE.time{fromsel}(FT_SPIKE.timestamp{fromsel}==curr(jj,2));
    
        toy=cidyymap(curr(jj+1,1));
        tosel=strcmp(FT_SPIKE.label,num2str(curr(jj+1,1)));
        tox=FT_SPIKE.time{tosel}(FT_SPIKE.timestamp{tosel}==curr(jj+1,2));
        
%         if ismember(in_maximum_tags{ii,1},{'SSL','BSL'})  % loop        
            
%         else
%             plot([fromx,tox],[fromy,toy],'-')
%         end
        plotmat=[plotmat;fromx,tox,fromy,toy];
    end

    plot(plotmat(:,1:2).',plotmat(:,3:4).','-','Color',cc(randi(size(cc,1)),:));
end
set(gca(),'XTick',5:0.2:6.2,'XTickLabel',0:0.2:1.2)
ylim([0.5,25.5])





