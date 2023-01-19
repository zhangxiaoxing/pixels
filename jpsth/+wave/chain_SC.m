if false
    load('chain_tag.mat','out');
end

memsess=-1;
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
greys=ephys.getGreyRegs('range','grey');

for dur=reshape(fieldnames(out),1,[])
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for cnid=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            % length selection
            if size(out.(dur{1}).(wv{1}).(cnid{1}).ts,2)<5
                continue
            end

            % region selection
            currsess=str2double(regexp(cnid,'(?<=s)[0-9]*(?=c)','match','once'));
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

            ccg_qual=out.(dur{1}).(wv{1}).(cnid{1}).ccg_qual;
            % 1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising
            % edge 6:falling edge
            strict_sel=~(ccg_qual(:,2)>=252) + ~(ccg_qual(:,4)>=2).*2 + ~(ccg_qual(:,4)<=40).*4; % + ~(ccg_qual(:,5)>245).*8;
            if all(strict_sel==0) % plot
                % TODO trials with multiple chains
                cntsid=out.(dur{1}).(wv{1}).(cnid{1}).ts_id;
                tsid1=cntsid(cntsid(:,3)==1,:);
                [~,tsidsel]=ismember(out.(dur{1}).(wv{1}).(cnid{1}).ts(:,1),tsid1(:,1));
                trls=tsid1(tsidsel,5);
                [gc,gr]=groupcounts(trls);

                % TODO plot raster
                for tt=reshape(find(gc>1),1,[])
                    figure('Position',[32,32,1440,240])
                    hold on;
                    for jj=1:size(out.(dur{1}).(wv{1}).(cnid{1}).ts,2)
                        ts=cntsid(cntsid(:,3)==jj & cntsid(:,5)==gr(tt),4)-1;
                        plot(ts,jj*ones(size(ts)),'k|')
                        % overlap chain activity
                        for kk=reshape(find(trls==gr(tt)),1,[])
                            suts=out.(dur{1}).(wv{1}).(cnid{1}).ts(kk,jj);
                            [~,tickpos]=ismember(suts,cntsid(cntsid(:,3)==jj & cntsid(:,5)==gr(tt),1));
                            plot(ts(tickpos),jj,'r|','LineWidth',2);
                        end
                    end
                    ylim([0.5,jj+0.5])
                    xlabel('Time (s)')
                    ylabel('SU #')
                    title("Dur "+dur{1}+", wave "+wv{1}+", #"+cnid{1}+", trial "+num2str(tt));
                    set(gca(),"YTick",1:jj,"YTickLabel",cnreg,'YDir','reverse')
                end

                % TODO plot ccgs
                if ~any(gc>1)
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

