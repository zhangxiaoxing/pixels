if false
    chainStr=load('chain_tag.mat','out');
elseif false
    bumpStr=load('chain_sust_tag_600.mat','out');
    out=bumpStr.out;
    chain_len_thres=6;
elseif true
    load('rings_wave_burst_600.mat','out')
    chain_len_thres=3;
end

skipccg=true;
warning("Skip ccg on matebook due to missing toolbox");

global_init;

memsess=-1;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
greys=ephys.getGreyRegs('range','grey');

for dur=reshape(fieldnames(out),1,[])
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for cnid=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            % length selection
            chain_len=numel(out.(dur{1}).(wv{1}).(cnid{1}).meta{1});
            if chain_len <chain_len_thres
                continue
            end

            % region selection
            if contains(cnid,'r')
                currsess=str2double(regexp(cnid,'(?<=s)[0-9]*(?=r)','match','once'));
            else
                currsess=str2double(regexp(cnid,'(?<=s)[0-9]*(?=c)','match','once'));
            end

            if currsess~=memsess
                sesssel=su_meta.sess==currsess;
                sesscid=su_meta.allcid(sesssel);
                sessreg=su_meta.reg_tree(5,sesssel).';
                memsess=currsess;
            end
            cncids=out.(dur{1}).(wv{1}).(cnid{1}).meta{1};
            [~,supos]=ismember(cncids,sesscid);
            cnreg=sessreg(supos);
            if ~all(ismember(cnreg,greys),'all') || numel(unique(cnreg))<2
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
                cntsid=out.(dur{1}).(wv{1}).(cnid{1}).ts_id;
                tsidOnset=arrayfun(@(x) cntsid(cntsid(:,3)==x,:),1:numel(cncids),'UniformOutput',false); 
                if isnumeric(out.(dur{1}).(wv{1}).(cnid{1}).ts) % one spike chain
                    [~,tsidsel]=ismember(out.(dur{1}).(wv{1}).(cnid{1}).ts(:,1),tsidOnset{1}(:,1));
                    [gc,gr]=groupcounts(trls);
                    gridx=reshape(find(gc>1),1,[]);
                    ttlist=gr(gridx);
                else % bump-wave
                    cndur=cellfun(@(x) diff(x([1,end],3),1,1),out.(dur{1}).(wv{1}).(cnid{1}).ts);
                    [maxdur,maxidx]=max(cndur);
                    if maxdur>7500 % 500 ms
                        onsetSUIdx=out.(dur{1}).(wv{1}).(cnid{1}).ts{maxidx}(1,1); % first SU index in chain/loops
                        pos_within=out.(dur{1}).(wv{1}).(cnid{1}).ts{maxidx}(1,1);
                        tsidsel=out.(dur{1}).(wv{1}).(cnid{1}).ts{maxidx}(1,2);
                        ttlist=tsidOnset{pos_within}(tsidsel,5); % current trial
                        % TODO: same trial, longest for each diff onset
                        % TODO: match trial for different 

                        cn_tsid_sel=cellfun(@(x) x(1,1:2),out.(dur{1}).(wv{1}).(cnid{1}).ts,'UniformOutput',false);
                        cn_trl_list=cellfun(@(x) tsidOnset{x(1)}(x(2),5),cn_tsid_sel);
                        onset_set=cell2mat(cn_tsid_sel(cn_trl_list==ttlist(1)).'); % [supos,pos_within_SU]

                        supp_onset=setdiff(unique(),tsidsel);
                        suppIdx=[];
                        for ii=1:numel(supp_onset)
                            onset_sel=find(cn_tsid_sel==supp_onset(ii));
                            [~,oneid]=max(cndur(onset_sel));
                            suppIdx=[suppIdx,onset_sel(oneid)];
                        end
                    else
                        tsidsel=[];
                        ttlist=[];
                    end
                end
                trls=tsidOnset(tsidsel,5);

                % TODO plot raster
                for tt=ttlist
                    figure('Position',[32,32,1440,240])
                    hold on;
                    for jj=1:chain_len
                        ts=cntsid(cntsid(:,3)==jj & cntsid(:,5)==tt,4)-1;
                        plot(ts,jj*ones(size(ts)),'|','Color',['#',dec2hex(jj,6)])
                        % overlap chain activity
                        if isnumeric(out.(dur{1}).(wv{1}).(cnid{1}).ts)
                            for kk=reshape(find(trls==tt),1,[])
                                suts=out.(dur{1}).(wv{1}).(cnid{1}).ts(kk,jj);
                                [~,tickpos]=ismember(suts,cntsid(cntsid(:,3)==jj & cntsid(:,5)==tt,1));
                                plot(ts(tickpos),jj,'|','LineWidth',2,'Color',['#FF',dec2hex(jj,4)]);
                            end
                        else
                            for iidx=[suppIdx,maxidx]
                                bts=out.(dur{1}).(wv{1}).(cnid{1}).ts{iidx};
                                % TODO: merge same trial different onset chains
                                suts=bts(bts(:,1)==jj,3);
                                [~,tickpos]=ismember(suts,cntsid(cntsid(:,3)==jj & cntsid(:,5)==tt,1));
                                plot(ts(tickpos),repmat(jj,numel(tickpos),1),'|','LineWidth',2,'Color',['#FF',dec2hex(jj,4)]);
                            end
                        end
                    end
                    ylim([0.5,jj+0.5])
                    xlabel('Time (s)')
                    ylabel('SU #')
                    title({"Dur "+dur{1}+", wave "+wv{1}+", #"+cnid{1}+", trial "+num2str(tt),...
                        num2str(cncids)});
                    set(gca(),"YTick",1:jj,"YTickLabel",cnreg,'YDir','reverse')
                    keyboard()
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
    end
end

