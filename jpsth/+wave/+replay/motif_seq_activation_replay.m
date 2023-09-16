function motif_seq_activation_replay(chain_replay,ring_replay,trials_dict)
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
sps=30000;

usess=union(chain_replay.session,ring_replay.session);
for sessid=reshape(usess,1,[])
    trials=cell2mat(trials_dict(sessid));
    before_sec=trials(1,1)./sps-60;
    session_tick=wave.replay.sessid2length(sessid);
    after_sec=(session_tick-trials(end,2))./sps-60-1;

    for wid=["s1d3","s2d3","s1d6","s2d6"]
        samp=str2double(regexp(wid,"(?<=s)\d(?=d)",'match','once')).*4;
        delay=str2double(regexp(wid,"(?<=d)\d$",'match','once'));
        trial_sel=setdiff(find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2)),1:5);
        ringsel=find(ring_replay.session==sessid & ring_replay.wave==wid);
        chainsel=find(chain_replay.session==sessid & chain_replay.wave==wid);
        for tt=reshape(trial_sel,1,[])
            trl_sel=intersect(find(trials(:,5)==samp),tt-5:tt);
            if numel(trl_sel)<3
                continue
            end
            pref_delay_sec=sum(diff(trials(trl_sel,1:2),1,2)./sps-1);

            np_trl_sel=intersect(find(trials(:,5)~=samp),tt-5:tt);
            np_delay_sec=sum(diff(trials(np_trl_sel,1:2),1,2)./sps-1);
            trials(end+1,:)=min(session_tick-3,trials(end,2)+14*sps);
            iti_sec=sum((trials(tt-4:tt+1,1)-trials(tt-5:tt,2))./sps-4); % 1s test + 3s response
            trials(end,:)=[];


            %% trial
            tonset=trials(tt,1);
            motif_spk=cell(0,1);
            motif_id=[];
            prefspk=0;
            npspk=0;
            itispk=0;
            beforespk=0;
            afterspk=0;
            for ridx=reshape(ringsel,1,[])
                trl_align=ring_replay.trl_align{ridx};
                plottrl=trl_align(:,1)>=tt-5 & trl_align(:,1)<=tt;% & all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                currts=cell2mat(ring_replay.ts{ridx}(plottrl));
                if ~isempty(currts)
                    motif_spk=[motif_spk;reshape(currts,1,[])];
                    motif_id=[motif_id;"r"+ridx];
                end
                prefspk=prefspk+nnz(ismember(trl_align(:,1),trl_sel) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1));
                npspk=npspk+nnz(ismember(trl_align(:,1),np_trl_sel) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1));
                itispk=itispk+nnz(ismember(trl_align(:,1),tt-5:tt) & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14)));
            end
            for cidx=reshape(chainsel,1,[])
                trl_align=chain_replay.trl_align{cidx};
                plottrl=trl_align(:,1)>=tt-5 & trl_align(:,1)<=tt;% & all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                currts=chain_replay.ts{cidx}(plottrl,:);
                if ~isempty(currts)
                    motif_spk=[motif_spk;reshape(currts,1,[])];
                    motif_id=[motif_id;"c"+cidx];
                end
                prefspk=prefspk+nnz(ismember(trl_align(:,1),trl_sel) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1));
                npspk=npspk+nnz(ismember(trl_align(:,1),np_trl_sel) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1));
                itispk=itispk+nnz(ismember(trl_align(:,1),tt-5:tt) & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14)));
            end

            if isempty(motif_id) || ~any(startsWith(motif_id,"r")) || ~any(startsWith(motif_id,"c"))
                continue
            end

            if numel(motif_id)<20 || numel(motif_id)>60
                continue
            end

            if ((prefspk./pref_delay_sec)/(npspk./np_delay_sec)<1.6)...
                    || ((prefspk./pref_delay_sec)/(itispk./iti_sec)<1.6)
                continue
            end
            
            
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
                    beforespk=beforespk+nnz(before_task);
                    currts=cell2mat(ring_replay.ts{ridx}(before_task));
                    if isempty(currts), currts={[]};end
                    before_motif_spk=[before_motif_spk;reshape(currts,1,[])];
                else % "c"
                    cidx=str2double(replace(mid,"c",""));
                    trl_align=chain_replay.trl_align{cidx};
                    if any(trl_align(:,9)>300),beforepass=true;end
                    before_task=trl_align(:,8)==1 & trl_align(:,9)<=120;
                    beforespk=beforespk+nnz(before_task);
                    currts=chain_replay.ts{cidx}(before_task,:);
                    if isempty(currts), currts={[]};end
                    before_motif_spk=[before_motif_spk;reshape(currts,1,[])];
                end
            end
            
            if ~beforepass || ((prefspk./pref_delay_sec)/(beforespk/before_sec)<1.6)
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
                    afterspk=afterspk+nnz(after_task);
                    currts=cell2mat(ring_replay.ts{ridx}(after_task));
                    if isempty(currts), currts={[]};end
                    after_motif_spk=[after_motif_spk;reshape(currts,1,[])];
                else % "c"
                    cidx=str2double(replace(mid,"c",""));
                    trl_align=chain_replay.trl_align{cidx};
                    if any(trl_align(:,2)>360),afterpass=true;end
                    after_task=trl_align(:,1)==size(trials,1) & trl_align(:,2)>60 & trl_align(:,2)<=180;
                    afterspk=afterspk+nnz(after_task);
                    currts=chain_replay.ts{cidx}(after_task,:);
                    if isempty(currts), currts={[]};end
                    after_motif_spk=[after_motif_spk;reshape(currts,1,[])];
                end
            end

            if ~afterpass || ((prefspk./pref_delay_sec)/(afterspk/after_sec)<1.6)
                continue
            end
            
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
            keyboard()
            appendfig("close",true,"multi",true,"fn","nested_loops_in_out_task.pdf","tag","in out task population motif activity")
            % close all
        end
    end
end



