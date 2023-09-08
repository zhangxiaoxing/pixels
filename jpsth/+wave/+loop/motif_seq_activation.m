function motif_seq_activation(chain_replay,ring_replay,trials_dict)
arguments
    chain_replay = [];
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
for sessid=36%reshape(usess,1,[])
    trials=cell2mat(trials_dict(sessid));
    for wid=["s1d3","s2d3","s1d6","s2d6"]
        samp=str2double(regexp(wid,"(?<=s)\d(?=d)",'match','once')).*4;
        delay=str2double(regexp(wid,"(?<=d)\d$",'match','once'));
        trial_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
        ringsel=find(ring_replay.session==sessid & ring_replay.wave==wid);
        chainsel=find(chain_replay.session==sessid & chain_replay.wave==wid);
        for tt=153%reshape(trial_sel,1,[])
            tonset=trials(tt,1);
            motif_spk=cell(0,1);
            motif_id=[];
            for ridx=reshape(ringsel,1,[])
                trl_align=ring_replay.trl_align{ridx};
                pref_delay=trl_align(:,1)==tt & all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                currts=cell2mat(ring_replay.ts{ridx}(pref_delay));
                if ~isempty(currts)
                    motif_spk=[motif_spk;reshape(currts,1,[])];
                    motif_id=[motif_id;"r"+ridx];
                end
            end
            for cidx=reshape(chainsel,1,[])
                trl_align=chain_replay.trl_align{cidx};
                pref_delay=trl_align(:,1)==tt & all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                currts=chain_replay.ts{cidx}(pref_delay,:);
                if ~isempty(currts)
                    motif_spk=[motif_spk;reshape(currts,1,[])];
                    motif_id=[motif_id;"c"+cidx];
                end
            end

            if numel(motif_id)<50 || numel(motif_id)>100
                continue
            end
            %% plot
            covered=false(1,3000);
            for rridx=1:size(motif_spk,1)
                raster=unique(round((motif_spk{rridx}-tonset-30000)./30)); % 30 @ 1ms resolution
                for ii=2:numel(raster)
                    if raster(ii)-raster(ii-1)<10
                        onset=ceil(raster(ii-1));
                        offset=ceil(raster(ii));
                        covered((onset+1):(offset+1))=true;
                    end
                end
            end
            
            if mean(covered)<0.20
                continue
            end

            edges = find(diff([0,covered,0]==1));
            eonset = edges(1:2:end-1);  % Start indices
            run_length = edges(2:2:end)-eonset;  % Consecutive ones counts

            if ~any(run_length>100)
                continue
            end

            allspk=cell2mat(motif_spk.');
            allspk_sec=(allspk-tonset-30000)./30000;
            sumhist=histcounts(allspk_sec,0:0.25:delay);

            per_motif_hist=cell2mat(cellfun(@(x) histcounts((x-tonset-30000)./30000,0:0.25:delay),motif_spk,'UniformOutput',false));
            TCOM=sum(repmat(0.125:0.25:delay,size(per_motif_hist,1),1).*per_motif_hist,2)./sum(per_motif_hist,2);
            [~,sidx]=sort(TCOM);
            cmap=colormap('lines');
            close all
            figure()
            tiledlayout(2,1);
            nexttile()
            hold on
            for rridx=1:size(motif_spk,1)
                raster=unique(round((motif_spk{sidx==rridx}-tonset-30000)./60)).*2; % 60 @ 2ms resolution
                scatter(raster./1000,rridx,'|','MarkerEdgeColor',cmap(rem(rridx,7)+1,:))
            end
            ymax=numel(motif_spk)+1;
            ylim([0,ymax])
          
            for ii=1:numel(eonset)
                % disp([eonset(ii),eonset(ii)+run_length(ii)]./1000)
                plot([eonset(ii),eonset(ii)+run_length(ii)]./1000,[ymax,ymax]-0.5,'k-','LineWidth',4);
            end
            xlim([0,delay])
            xlabel("Time (s)")
            ylabel("Motif #")
            title("S"+sessid+"T"+tt)

            nexttile()
            plot(0.125:0.25:delay,sumhist,'r-','LineWidth',1);
            xlabel("Time (s)")
            ylabel("Motif spike frequency (Hz)")
            appendfig("close",true,"fn","per_trl_motif.pdf");
        end
    end
end



