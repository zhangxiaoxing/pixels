statstype="chain";
skipccg=false;
skipraster=false;
chain_len_thres=4;
accu_spk_thres=3;

if ~exist('chain_replay','var')
switch statstype
    case "chain"
        memstr=load(fullfile("binary","motif_replay.mat"),'chain_replay');
        nmstr=load(fullfile('binary','motif_replay_chain_nonmem.mat'),'chain_replay');
        chain_replay=[memstr.chain_replay;nmstr.chain_replay];

        load(fullfile('binary','trials_dict.mat'),'trials_dict');
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

for cid=[2405,2631]
% for cid=1:size(chain_replay,1)
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
        strict_sel=~(ccg_qual(:,2)>=252 & ccg_qual(:,2)<=265) + ~(ccg_qual(:,4)>=3).*2 + ~(ccg_qual(:,4)<=30).*4; % + ~(ccg_qual(:,5)>245).*8;
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
        sgtitle("S"+chain_replay.session(cid)+"W"+wid+"C"+cid+"N"+max(gc));
        % appendfig('tag',['chain showcase',strjoin(cnreg.')])
        if skipraster
            keyboard()
        else
            [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(chain_replay.session(cid),'keep_trial',true);
            plot_trl_sel=gr(gc>=accu_spk_thres);

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
                title({"S"+chain_replay.session(cid)+"C"+cid;strjoin(cnreg,'>')})
            end
            % keyboard()
        end
    end
end
appendfig('fn','chain_SC.pdf','tag','Chain showcase','multi',true)