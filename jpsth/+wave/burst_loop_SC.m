%specifically burst loop
accu_spk_thres=16;
loop_len_thres=4;


dbfile=fullfile("bzdata","rings_wave_burst_iter_150.db");
conn=sqlite(dbfile,'readonly');
keys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
% conn.close();

tskeys=keys(endsWith(keys,"_ts"));

skipccg=true;
warning("Skip ccg on matebook due to missing toolbox");

global_init;
memsess=-1;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
greys=ephys.getGreyRegs('range','grey');

% for dur=reshape(fieldnames(out),1,[])
%     for wv=reshape(fieldnames(out.(dur{1})),1,[])
%         for cnid=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
% length selection

for tskey=reshape(tskeys,1,[])
    metakey=replace(tskey,"_ts","_meta");
    loopcids=table2array(conn.sqlread(metakey));
    loop_len=numel(loopcids);
    if loop_len < loop_len_thres
        continue
    end
    loop_ts=table2array(conn.sqlread(tskey));
    [accu_spk,occ_id]=groupcounts(loop_ts(:,1));
    if max(accu_spk)<accu_spk_thres % split checks due to performance concern
        continue
    end

    % region selection
    currsess=str2double(regexp(tskey,'(?<=s)\d+(?=r)','match','once'));
    if currsess~=memsess
        sesssel=su_meta.sess==currsess;
        sesscid=su_meta.allcid(sesssel);
        sessreg=su_meta.reg_tree(5,sesssel).';
        memsess=currsess;
    end

    [~,supos]=ismember(loopcids,sesscid);
    cnreg=sessreg(supos);
    if ~all(ismember(cnreg,greys),'all') || numel(unique(cnreg))<2 %TODO: within-region?
        continue
    end
    if skipccg
        strict_sel=0;
    else
        ccg_qual=out.(dur{1}).(wv{1}).(cnid{1}).ccg_qual;
        % 1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising
        % edge 6:falling edge
        strict_sel=~(ccg_qual(:,2)>=252) + ~(ccg_qual(:,4)>=2).*2 + ~(ccg_qual(:,4)<=40).*4; % + ~(ccg_qual(:,5)>245).*8;
    end
    if all(strict_sel==0) % plot
        loop_tsid=table2array(conn.sqlread(replace(tskey,"_ts","_tsid")));
        ts_id_per_cid=arrayfun(@(x) loop_tsid(loop_tsid(:,3)==x,:),1:numel(loopcids), ...
            'UniformOutput',false);
        tsidOnset=splitapply(@(x) x(1,:), loop_ts(:,2:4),loop_ts(:,1));
        
        loop_trl_list=arrayfun(@(x) ts_id_per_cid{tsidOnset(x,1)}(tsidOnset(x,2),5),1:size(tsidOnset,1));% corresponding trial

        trl_su_tspos=unique(cell2mat(arrayfun(@(x) [loop_trl_list(x),tsidOnset(x,:)],(1:size(tsidOnset,1)).','UniformOutput',false)),'rows'); % [trial, in-chain-idx, per-su-idx]

        tsdiff=diff(trl_su_tspos,1,1);

        trl_su_tspos(tsdiff(:,1)==0 & tsdiff(:,2)==0 & tsdiff(:,4)<150,:)=[];

        [gc,gr]=groupcounts(trl_su_tspos(:,1));


        % plot raster
        for tt=(gr(gc>1)).'
            trl_ts=find(loop_trl_list==tt);
            join_ts=loop_ts(ismember(loop_ts(:,1),trl_ts),2:4);
            figure('Position',[32,32,1440,240])
            hold on;
            for jj=1:loop_len
                ts=loop_tsid(loop_tsid(:,3)==jj & loop_tsid(:,5)==tt,4)-1;
                plot(ts,jj*ones(size(ts)),'|','Color',['#',dec2hex(jj,6)])
                % overlap chain activity
                suts=unique(join_ts(join_ts(:,1)==jj,3));
                [~,tickpos]=ismember(suts,loop_tsid(loop_tsid(:,3)==jj & loop_tsid(:,5)==tt,1));
                plot(ts(tickpos),repmat(jj,numel(tickpos),1),'|','LineWidth',2,'Color',['#FF',dec2hex(jj,4)]);
            end
            ylim([0.5,jj+0.5])
            xlabel('Time (s)')
            ylabel('SU #')
            title({replace(tskey,"_ts","")+", trial "+num2str(tt),...
                num2str(loopcids)});
            set(gca(),"YTick",1:jj,"YTickLabel",cnreg,'YDir','reverse')
            gh=groot;
            if numel(gh.Children)>=20
                keyboard()
            end
        end

        % TODO plot ccgs
        if skipccg || ~any(gc>1)
            continue
        end
        figure()
        tiledlayout('flow')
        for ii=1:size(out.(dur{1}).(wv{1}).(cnid{1}).ccgs,1)
            nexttile()
            hold on
            plot(out.(dur{1}).(wv{1}).(cnid{1}).ccgs(ii,201:301))
            xline(51,'--k')
            xlim([1,101]);
        end
        sgtitle("Dur "+dur{1}+", wave "+wv{1}+", #"+cnid{1});
        keyboard();

    else
        disp(strict_sel)
        %                 keyboard()
    end
end
