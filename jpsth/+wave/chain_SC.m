statstype="chain";
skipccg=false;
if false
switch statstype
    case "chain"
        load(fullfile('binary','motif_replay.mat'),'chain_replay');
        load(fullfile('binary','trials_dict.mat'),'trials_dict');
        chain_len_thres=4;
        accu_spk_thres=3;
    case"burstchain"
        fstr=load(fullfile('bzdata','chain_sust_tag_150.mat'),'out');
        out=fstr.out;
        chain_len_thres=5;
        accu_spk_thres=10;
    otherwise
        disp("data type mismatch")
        keyboard()
end
sps=30000;
global_init;

memsess=-1;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
load(fullfile("binary","sums_conn_10.mat"),'sums_conn_str');
greys=ephys.getGreyRegs('range','grey');

end
for cid=2388 %1:size(chain_replay,1)
    if ~chain_replay.meta{cid,3}
        continue
    end
    cncids=chain_replay.meta{cid,2};
    chain_len=numel(cncids);
    if chain_len <chain_len_thres
        continue
    end

    if statstype=="chain"
        wid=chain_replay.wave(cid);
        % samp=str2double(regexp(wid,"(?<=s)\d(?=d)",'match','once')).*4;
        delay=str2double(regexp(wid,"(?<=d)\d$",'match','once'));
        trials=cell2mat(trials_dict(chain_replay.session(cid)));
        % trl_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
        % pref_delay_sec=sum(diff(trials(trl_sel,1:2),1,2)./sps-1);
        
        trl_align=chain_replay.trl_align{cid};
        pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
        if ~any(pref_delay)
            continue
        end
        [gc,gr]=groupcounts(trl_align(pref_delay,1));

        if max(gc)<accu_spk_thres
            continue
        end
        % TODO: skipped for now
        % maxcnts=out.(dur{1}).(wv{1}).(cnid{1}).ts{max_ts};
    end

    % region selection
    if chain_replay.session(cid)~=memsess
        sesssel=su_meta.sess==chain_replay.session(cid);
        sesscid=su_meta.allcid(sesssel);
        sessreg=su_meta.reg_tree(5,sesssel).';
        memsess=chain_replay.session(cid);
        [spkID,~,~,~,~,~]=ephys.getSPKID_TS(chain_replay.session(cid),'keep_trial',false);
    end

    [~,supos]=ismember(cncids,sesscid);
    cnreg=sessreg(supos);

    if skipccg
        strict_sel=0;
    else
        %TODO ccg
        spath=replace(ephys.sessid2path(chain_replay.session(cid)),'\','/');
        ccgssel=find(contains({sums_conn_str.folder},spath));
        sigcid=sums_conn_str(ccgssel).sig_con;
        edges=[cncids(1:end-1);cncids(2:end)].';
        ccg_qual=[];
        ccg_seq=[];
        for ii=1:size(edges,1)
            ccgsel=find(sigcid(:,1)==edges(ii,1) & sigcid(:,2)==edges(ii,2));
            ccg_qual=[ccg_qual;sums_conn_str(ccgssel).qc(ccgsel,:)];
            ccg_seq=[ccg_seq,ccgsel];
        end
        % 1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising
        % edge 6:falling edge
        strict_sel=~(ccg_qual(:,2)>=253 & ccg_qual(:,2)<=265) + ~(ccg_qual(:,4)>=2).*2 + ~(ccg_qual(:,4)<=35).*4; % + ~(ccg_qual(:,5)>245).*8;
    end
    if all(strict_sel==0) % plot
        %%
        figure()
        tiledlayout('flow')
        for ii=1:numel(ccg_seq)
            nexttile()
            hold on
            normfactor=2500./nnz(spkID==cncids(ii)); % Hz in 0.4ms bin, 1000/0.4
            % normfactor=1;
            plot(-100:0.4:100,sums_conn_str(ccgssel).ccg_sc(ccg_seq(ii),:).*normfactor,'-r')
            xline([-10,0,10,],'--k')
            xlim([-12,24]);
            xlabel('Latency (ms)')
            ylabel('Firing rate (Hz)')
        end
        sgtitle("s"+chain_replay.session(cid)+"w"+wid+"C"+cid);
        % appendfig('tag',['chain showcase',strjoin(cnreg.')])
        keyboard()
        [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(chain_replay.session(cid),'keep_trial',true);
        plot_trl_sel=gr(gc>4);
        
        for currtrl=reshape(plot_trl_sel,1,{})
            plot_sel=trl_align(:,1)==currtrl & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            trlonset=trials(currtrl,1);
            tsall=(chain_replay.ts{cid}(plot_sel,:)-trlonset-30000)./30000;


            figure()
            cmap=colormap('lines');
            hold on;
            for pp=1:size(tsall,2)
                ftsel=strcmp(FT_SPIKE.label,num2str(cncids(pp)));
                otherspk=FT_SPIKE.time{ftsel}(FT_SPIKE.trial{ftsel}==currtrl & FT_SPIKE.time{ftsel}>=1 & FT_SPIKE.time{ftsel}<delay+1);
                plot([otherspk;otherspk],repmat([pp-0.4,pp+0.4],numel(otherspk),1).','-','Color',[0.5,0.5,0.48+pp/255],'LineWidth',0.5)
                plot([tsall(:,pp),tsall(:,pp)].',repmat([pp-0.4,pp+0.4],size(tsall,1),1).','-','Color',cmap(pp,:),'LineWidth',1.5)
            end
            set(gca,'YDir','reverse','XTick',1:2:7,'XTickLabel',0:2:6)
            ylim([0.5,numel(cncids)+0.5])
        end
    end
end

function others()
        %%
        for ii=[]
            
        cn_tsid=out.(dur{1}).(wv{1}).(cnid{1}).ts_id;
        tsid_per_cid=arrayfun(@(x) cn_tsid(cn_tsid(:,3)==x,:),1:numel(cncids),'UniformOutput',false);
        if isnumeric(out.(dur{1}).(wv{1}).(cnid{1}).ts) % one spike chain
            [~,tsidsel]=ismember(out.(dur{1}).(wv{1}).(cnid{1}).ts(:,1),tsid_per_cid{1}(:,1));
            cn_trl_list=tsid_per_cid{1}(tsidsel,5);
            [gc,gr]=groupcounts(cn_trl_list);
            %                     gridx=reshape(find(gc>1),1,[]);
            %                     ttlist=gr(gridx);
        else % bursting wave

            cn_first_tagged=cellfun(@(x) x(1,1:3),out.(dur{1}).(wv{1}).(cnid{1}).ts,'UniformOutput',false); % first tagged spike, in chain idx and per-su-ts-idx
            cn_trl_list=cellfun(@(x) tsid_per_cid{x(1)}(x(2),5),cn_first_tagged);% corresponding trial
            trl_su_tspos=unique(cell2mat(arrayfun(@(x) [cn_trl_list(x),cn_first_tagged{x}],(1:numel(cn_trl_list)).','UniformOutput',false)),'rows'); % [trial, in-chain-idx, per-su-idx]

            burst_trl=cn_trl_list(per_seq_spk>=accu_spk_thres);
            trl_su_tspos=trl_su_tspos(ismember(trl_su_tspos(:,1),burst_trl),:);

            tsdiff=diff(trl_su_tspos,1,1);
            trl_su_tspos(tsdiff(:,1)==0 & tsdiff(:,2)==0 & tsdiff(:,3)==1 & tsdiff(:,4)<300,:)=[];

            [gc,gr]=groupcounts(trl_su_tspos(:,1));
            %                     maxtrl=cn_trl_list(max_ts);
            %                     if gc(gr==maxtrl)<2
            %                         continue
            %                     end
        end

        % plot raster
        for tt=(gr(gc>1)).'
            if statstype~="chain"
                trl_ts=cellfun(@(x) tsid_per_cid{x(1)}(x(2),5),cn_first_tagged)==tt;
                join_ts=cell2mat(out.(dur{1}).(wv{1}).(cnid{1}).ts(trl_ts).');
            end

            figure('Position',[32,32,1440,240])
            hold on;
            for jj=1:chain_len
                ts=cn_tsid(cn_tsid(:,3)==jj & cn_tsid(:,5)==tt,4)-1;
                plot(ts,jj*ones(size(ts)),'|','Color',['#',dec2hex(jj,6)])
                % overlap chain activity
                if isnumeric(out.(dur{1}).(wv{1}).(cnid{1}).ts)
                    for kk=reshape(find(cn_trl_list==tt),1,[])
                        suts=out.(dur{1}).(wv{1}).(cnid{1}).ts(kk,jj);
                        [~,tickpos]=ismember(suts,cn_tsid(cn_tsid(:,3)==jj & cn_tsid(:,5)==tt,1));
                        plot(ts(tickpos),jj,'|','LineWidth',2,'Color',['#FF',dec2hex(jj,4)]);
                    end
                else
                    suts=unique(join_ts(join_ts(:,1)==jj,3));
                    [~,tickpos]=ismember(suts,cn_tsid(cn_tsid(:,3)==jj & cn_tsid(:,5)==tt,1));
                    plot(ts(tickpos),repmat(jj,numel(tickpos),1),'|','LineWidth',2,'Color',['#FF',dec2hex(jj,4)]);

                    %                             [~,tickpos]=ismember(maxcnts(maxcnts(:,1)==jj,3),cn_tsid(cn_tsid(:,3)==jj & cn_tsid(:,5)==tt,1));
                    %                             plot(ts(tickpos),repmat(jj,numel(tickpos),1),'|','LineWidth',2,'Color',['#',dec2hex(jj,2),'FF00']);

                end
            end
            ylim([0.5,jj+0.5])
            xlabel('Time (s)')
            ylabel('SU #')
            title({"Dur "+dur{1}+", wave "+wv{1}+", #"+cnid{1}+", trial "+num2str(tt),...
                num2str(cncids)});
            set(gca(),"YTick",1:jj,"YTickLabel",cnreg,'YDir','reverse')
            gh=groot;
            if numel(gh.Children)>=50
                keyboard()
            end
        end

        % TODO plot ccgs
        if skipccg || ~any(gc>1)
            continue
        end

        
        % CCG PLOT
        appendfig('fn','chain_sc.pdf','close',true,'multi',true,'path',fullfile('Z:'))
    % if
    % else
    %     %                 disp(strict_sel)
    %     %                 keyboard()
    end
end


function chain_replay_SC()
%%
% sschain=load(fullfile('bzdata','chain_tag_tbl.mat'),'out');
[sschain.out,unfound]=wave.chain_tag.tag(chains_uf,'skip_save',true,'len_thresh',len_thresh,'odor_only',false,'extend_trial',true,'sesses',100,'idces',1312); % per-spk association
out=sschain.out;
dur={'d3'};
wv={'olf_s1'};
cnid={'s100c1312'};
statstype="chain";

cn_tsid=out.(dur{1}).(wv{1}).(cnid{1}).ts_id;
cncids=out.(dur{1}).(wv{1}).(cnid{1}).meta{1};
tsid_per_cid=arrayfun(@(x) cn_tsid(cn_tsid.POS==x,:),1:numel(cncids),'UniformOutput',false);
[~,tsidsel]=ismember(out.(dur{1}).(wv{1}).(cnid{1}).ts(:,1),tsid_per_cid{1}.TS);
cn_trl_list=tsid_per_cid{1}.Trial(tsidsel);

su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
sesssel=su_meta.sess==100;
sesscid=su_meta.allcid(sesssel);
sessreg=su_meta.reg_tree(5,sesssel).';
[~,supos]=ismember(cncids,sesscid);
cnreg=sessreg(supos);

for tt=204
    if statstype~="chain"
        trl_ts=cellfun(@(x) tsid_per_cid{x(1)}(x(2),5),cn_first_tagged)==tt;
        join_ts=cell2mat(out.(dur{1}).(wv{1}).(cnid{1}).ts(trl_ts).');
    end

    figure('Position',[32,32,1440,240])
    hold on;
    for jj=1:numel(cncids)
        ts=cn_tsid.Time(cn_tsid.POS==jj & cn_tsid.Trial==tt)-1;
        plot(ts,jj*ones(size(ts)),'|','Color',['#',dec2hex(jj,6)])
        % overlap chain activity
        if isnumeric(out.(dur{1}).(wv{1}).(cnid{1}).ts)
            for kk=reshape(find(cn_trl_list==tt),1,[])
                suts=out.(dur{1}).(wv{1}).(cnid{1}).ts(kk,jj);
                [~,tickpos]=ismember(suts,cn_tsid.TS(cn_tsid.POS==jj & cn_tsid.Trial==tt));
                plot(ts(tickpos),jj,'|','LineWidth',2,'Color',['#FF',dec2hex(jj,4)]);
            end
        else
            suts=unique(join_ts(join_ts(:,1)==jj,3));
            [~,tickpos]=ismember(suts,cn_tsid.TS(cn_tsid.POS==jj & cn_tsid.Trial==tt));
            plot(ts(tickpos),repmat(jj,numel(tickpos),1),'|','LineWidth',2,'Color',['#FF',dec2hex(jj,4)]);
        end
    end
    ylim([0.5,jj+0.5])
    xlabel('Time (s)')
    ylabel('SU #')
    title({"Dur "+dur{1}+", wave "+wv{1}+", #"+cnid{1}+", trial "+num2str(tt),...
        num2str(cncids)});
    set(gca(),"YTick",1:jj,"YTickLabel",cnreg,'YDir','reverse')
    gh=groot;
    if numel(gh.Children)>=50
        keyboard()
    end
end
%%
end

