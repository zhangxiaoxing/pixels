statstype="chain";
skipccg=false;

switch statstype
    case "chain"
        load(fullfile('bzdata','chain_tag.mat'),'out');
        chain_len_thres=5;
        accu_spk_thres=5;
    case"burstchain"
        fstr=load(fullfile('bzdata','chain_sust_tag_150.mat'),'out');
        out=fstr.out;
        chain_len_thres=5;
        accu_spk_thres=10;
    otherwise
        disp("data type mismatch")
        keyboard()
end

if skipccg
    warning("Skip ccg on matebook due to missing toolbox");
end

global_init;

memsess=-1;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
greys=ephys.getGreyRegs('range','grey');

for dur=reshape(fieldnames(out),1,[])
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for cnid=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            % length selection
            cncids=out.(dur{1}).(wv{1}).(cnid{1}).meta{1};
            chain_len=numel(cncids);
            if chain_len <chain_len_thres
                continue
            end
            if statstype~="chain"
                per_seq_spk=cellfun(@(x) size(x,1),out.(dur{1}).(wv{1}).(cnid{1}).ts);
                [spkcount,max_ts]=max(per_seq_spk);
                if spkcount<accu_spk_thres
                    continue
                end
                maxcnts=out.(dur{1}).(wv{1}).(cnid{1}).ts{max_ts};
            end

            % region selection
            currsess=str2double(regexp(cnid,'(?<=s)[0-9]*(?=c)','match','once'));
            if currsess~=memsess
                sesssel=su_meta.sess==currsess;
                sesscid=su_meta.allcid(sesssel);
                sessreg=su_meta.reg_tree(5,sesssel).';
                memsess=currsess;
            end
            
            [~,supos]=ismember(cncids,sesscid);
            cnreg=sessreg(supos);
            if ~all(ismember(cnreg,greys),'all') % || numel(unique(cnreg))<2
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

                [spkID,~,~,~,~,~]=ephys.getSPKID_TS(currsess,'keep_trial',false);
                figure()
                tiledlayout('flow')
                for ii=1:size(out.(dur{1}).(wv{1}).(cnid{1}).ccgs,1)
                    nexttile()
                    hold on
                    normfactor=2500./nnz(spkID==out.(dur{1}).(wv{1}).(cnid{1}).meta{1}(ii)); % Hz in 0.4ms bin, 1000/0.4
                    plot(-100:0.4:100,out.(dur{1}).(wv{1}).(cnid{1}).ccgs(ii,:).*normfactor,'-r')
                    xline([-10,0,10,],'--k')
                    xlim([-12,24]);
                    xlabel('Latency (ms)')
                    ylabel('Firing rate (Hz)')
                end
                sgtitle("Dur "+dur{1}+", wave "+wv{1}+", #"+cnid{1});
                appendfig('fn','chain_sc.pdf','close',true,'multi',true,'path',fullfile('Z:'))

            else
%                 disp(strict_sel)
%                 keyboard()
            end
        end
    end
end

