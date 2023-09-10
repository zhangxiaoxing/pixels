function motif_seq_activation(chain_replay,ring_replay,trials_dict)
arguments
    chain_replay = []
    ring_replay = []
    trials_dict = []
end
if isempty(chain_replay) || isempty (ring_replay)
    load(fullfile('binary','motif_replay.mat'),'ring_replay','chain_replay');
end
if isempty(trials_dict)
    load(fullfile('binary','trials_dict.mat'),'trials_dict');
end

usess=union(chain_replay.session,ring_replay.session);
for sessid=reshape(usess,1,[])
    trials=cell2mat(trials_dict(sessid));
    for wid=["s1d3","s2d3","s1d6","s2d6"]
        samp=str2double(regexp(wid,"(?<=s)\d(?=d)",'match','once')).*4;
        delay=str2double(regexp(wid,"(?<=d)\d$",'match','once'));
        trial_sel=setdiff(find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2)),1:5);
        ringsel=find(ring_replay.session==sessid & ring_replay.wave==wid);
        chainsel=find(chain_replay.session==sessid & chain_replay.wave==wid);
        for tt=reshape(trial_sel,1,[])
            %% trial
            tonset=trials(tt,1);
            motif_spk=cell(0,1);
            motif_id=[];
            for ridx=reshape(ringsel,1,[])
                trl_align=ring_replay.trl_align{ridx};
                pref_delay=trl_align(:,1)>=tt-5 & trl_align(:,1)<=tt;% & all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                currts=cell2mat(ring_replay.ts{ridx}(pref_delay));
                if ~isempty(currts)
                    motif_spk=[motif_spk;reshape(currts,1,[])];
                    motif_id=[motif_id;"r"+ridx];
                end
            end
            for cidx=reshape(chainsel,1,[])
                trl_align=chain_replay.trl_align{cidx};
                pref_delay=trl_align(:,1)>=tt-5 & trl_align(:,1)<=tt;% & all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                currts=chain_replay.ts{cidx}(pref_delay,:);
                if ~isempty(currts)
                    motif_spk=[motif_spk;reshape(currts,1,[])];
                    motif_id=[motif_id;"c"+cidx];
                end
            end

            if isempty(motif_id) || ~any(startsWith(motif_id,"r")) || ~any(startsWith(motif_id,"c"))
                continue
            end

            if numel(motif_id)<40 || numel(motif_id)>100
                continue
            end

            % covered=false(1,3000);
            % for rridx=1:size(motif_spk,1)
            %     raster=unique(round((motif_spk{rridx}-tonset-30000)./30)); % 30 @ 1ms resolution
            %     for ii=2:numel(raster)
            %         if raster(ii)-raster(ii-1)<10
            %             onset=ceil(raster(ii-1));
            %             offset=ceil(raster(ii));
            %             covered((onset+1):(offset+1))=true;
            %         end
            %     end
            % end
            % 
            % if mean(covered)<0.20
            %     continue
            % end

            % edges = find(diff([0,covered,0]==1));
            % eonset = edges(1:2:end-1);  % Start indices
            % run_length = edges(2:2:end)-eonset;  % Consecutive ones counts
            % 
            % if ~any(run_length>100)
            %     continue
            % end

            % allspk=cell2mat(motif_spk.');
            % allspk_sec=(allspk-tonset-30000)./30000;
            % sumhist=histcounts(allspk_sec,0:0.25:delay);

            
            %% before 

            before_onset=trials(1,1);
            before_motif_spk=cell(0,1);
            beforepass=false;
            for mid=reshape(motif_id,1,[])
                if startsWith(mid,"r")
                    ridx=str2double(replace(mid,"r",""));
                    trl_align=ring_replay.trl_align{ridx};
                    if any(trl_align(:,9)>300),beforepass=true;end
                    before_task=trl_align(:,8)==1 & trl_align(:,9)<=120;
                    currts=cell2mat(ring_replay.ts{ridx}(before_task));
                    if isempty(currts), currts={[]};end
                    before_motif_spk=[before_motif_spk;reshape(currts,1,[])];
                else % "c"
                    cidx=str2double(replace(mid,"c",""));
                    trl_align=chain_replay.trl_align{cidx};
                    if any(trl_align(:,9)>300),beforepass=true;end
                    before_task=trl_align(:,8)==1 & trl_align(:,9)<=120;
                    currts=chain_replay.ts{cidx}(before_task,:);
                    if isempty(currts), currts={[]};end
                    before_motif_spk=[before_motif_spk;reshape(currts,1,[])];
                end
            end
            
            if ~beforepass
                continue
            end

            %% after 
            after_onset=trials(end,1);
            after_motif_spk=cell(0,1);
            afterpass=false;
            for mid=reshape(motif_id,1,[])
                if startsWith(mid,"r")
                    ridx=str2double(replace(mid,"r",""));
                    trl_align=ring_replay.trl_align{ridx};
                    if any(trl_align(:,2)>360),afterpass=true;end
                    after_task=trl_align(:,1)==size(trials,1) & trl_align(:,2)>60 & trl_align(:,2)<=180;
                    currts=cell2mat(ring_replay.ts{ridx}(after_task));
                    if isempty(currts), currts={[]};end
                    after_motif_spk=[after_motif_spk;reshape(currts,1,[])];
                else % "c"
                    cidx=str2double(replace(mid,"c",""));
                    trl_align=chain_replay.trl_align{cidx};
                    if any(trl_align(:,2)>360),afterpass=true;end
                    after_task=trl_align(:,1)==size(trials,1) & trl_align(:,2)>60 & trl_align(:,2)<=180;
                    currts=chain_replay.ts{cidx}(after_task,:);
                    if isempty(currts), currts={[]};end
                    after_motif_spk=[after_motif_spk;reshape(currts,1,[])];
                end
            end

            if ~afterpass
                continue
            end
            
            % selectivity criteria
            for shadet=tt-5:tt
                if trials(shadet,5)==trials(tt,5)
                else
                end
            end



            % before after half of pref

            % nonpref half of pref

            % ITI half of pref



            %% plot trial
            % skip sort due to low visibiliby
            % per_motif_hist=cell2mat(cellfun(@(x) histcounts((x-tonset-30000)./30000,0:0.25:delay),motif_spk,'UniformOutput',false));
            % TCOM=sum(repmat(0.125:0.25:delay,size(per_motif_hist,1),1).*per_motif_hist,2)./sum(per_motif_hist,2);
            % [~,sidx]=sort(TCOM);
            fh=figure(); % avoid use tiled layout for performance reasons.
            cmap=colormap('lines');
            hold on

            % shade
            for shadet=tt-5:tt
                fill([trials(shadet,1),trials(shadet,1)+30000,trials(shadet,1)+30000,trials(shadet,1)]./30000,...
                    [0,0,numel(motif_id)+1,numel(motif_id)+1],...
                    'k','EdgeColor','none','FaceAlpha',0.1)
                fill([trials(shadet,2),trials(shadet,2)+30000,trials(shadet,2)+30000,trials(shadet,2)]./30000,...
                    [0,0,numel(motif_id)+1,numel(motif_id)+1],...
                    'b','EdgeColor','none','FaceAlpha',0.1)
                if trials(shadet,5)==trials(tt,5)
                    fill([trials(shadet,1)+30000,trials(shadet,2),trials(shadet,2),trials(shadet,1)+30000]./30000,...
                    [0,0,numel(motif_id)+1,numel(motif_id)+1],...
                    'r','EdgeColor','none','FaceAlpha',0.1)
                else
                    fill([trials(shadet,1)+30000,trials(shadet,2),trials(shadet,2),trials(shadet,1)+30000]./30000,...
                    [0,0,numel(motif_id)+1,numel(motif_id)+1],...
                    'y','EdgeColor','none','FaceAlpha',0.1)
                end
            end

            % raster

            for rridx=1:size(motif_spk,1)
                raster=unique(round((motif_spk{rridx})./600)).*20; % 600 @ 20ms resolution
                scatter(raster./1000,rridx,'|','MarkerEdgeColor',cmap(rem(rridx,7)+1,:))
            end
            ymax=numel(motif_spk)+1;
            ylim([0,ymax])
          
            % for ii=1:numel(eonset)
            %     % disp([eonset(ii),eonset(ii)+run_length(ii)]./1000)
            %     plot([eonset(ii),eonset(ii)+run_length(ii)]./1000,[ymax,ymax]-0.5,'k-','LineWidth',4);
            % end
            xlim([trials(tt-5,1),trials(tt+1,1)]./30000)
            xlabel("Time (s)")
            ylabel("Motif #")
            title("S"+sessid+"T"+tt)


            %% plot before
            fh=figure();
            cmap=colormap('lines');
            hold on
            for rridx=1:size(before_motif_spk,1)
                disp(rridx);
                if isempty(before_motif_spk{rridx}),continue;end
                raster=unique(round((before_motif_spk{rridx}-before_onset)./600)).*20; % 600 @ 20ms resolution
                scatter(raster./1000,rridx,'|','MarkerEdgeColor',cmap(rem(rridx,7)+1,:))
            end
            ymax=numel(before_motif_spk)+1;
            ylim([0,ymax])
            xlim([-120,0])
            xlabel("Time (s)")
            ylabel("Motif #")
            title("S"+sessid+"T"+tt+"Before")

            %% plot after
            fh=figure();
            cmap=colormap('lines');
            hold on
            for rridx=1:size(after_motif_spk,1)
                disp(rridx);
                if isempty(after_motif_spk{rridx}),continue;end
                raster=unique(round((after_motif_spk{rridx}-after_onset)./600)).*20; % 1500 @ 50ms resolution
                scatter(raster./1000,rridx,'|','MarkerEdgeColor',cmap(rem(rridx,7)+1,:))
            end
            ymax=numel(after_motif_spk)+1;
            ylim([0,ymax])
            xlim([60,180])
            xlabel("Time (s)")
            ylabel("Motif #")
            title("S"+sessid+"T"+tt+"After")
            appendfig("close",true,"multi",true,"fn","nested_loops_in_out_task.pdf","tag","in out task population motif activity")
            close all
        end
    end
end



