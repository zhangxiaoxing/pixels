% TODO: rearrange the functions to facilitate the combination of 6s and 3s
% stats when beneficial

classdef fc_com_reg_wave < handle
    methods (Static)
        function fc_com_pvsst_stats=stats(wave_meta,su_com_map,opt)
            arguments
                wave_meta
                su_com_map
                opt.delay (1,1) double {mustBeMember(opt.delay,[3 6])} = 6
                opt.odor_only (1,1) logical = false
            end

            [sig,~]=bz.load_sig_sums_conn_file();
            sig=bz.join_fc_waveid(sig,wave_meta.wave_id);
            usess=unique(sig.sess);
            fc_com_pvsst_stats=[];

            for sii=reshape(usess,1,[]) %iterate through sessions
                sesssel=sig.sess==sii;

                suid=sig.suid(sesssel,:);
                waveid=sig.waveid(sesssel,:);
                regsess=squeeze(sig.reg(sesssel,5,:));

                com_sess_pct=nan(size(suid));
                if opt.odor_only
                    ffs=["olf_s1","olf_s2"];
                    durkey='com';
                else
                    if opt.delay==6
                        ffs=["s1d6","s2d6","olf_s1","olf_s2","dur_d6"];
                        durkey='com6';
                    elseif opt.delay==3
                        ffs=["s1d3","s2d3","olf_s1","olf_s2","dur_d3"];
                        durkey='com3';
                    end
                end
                for ff=ffs
                    if isfield(su_com_map.(['s',num2str(sii)]),ff)
                        sukeys=su_com_map.(['s',num2str(sii)]).(ff).(durkey).keys(); % prefered SUid
                        susel=ismember(suid,int32(cell2mat(sukeys)));% same dim as suid
                        com_sess_pct(susel)=cell2mat(su_com_map.(['s',num2str(sii)]).(ff).(durkey).values(num2cell(suid(susel)))); % out put is nx2 in dim
                    end
                end
                fc_com_pvsst_stats=[fc_com_pvsst_stats;double(sii).*ones(size(suid(:,1))),double(suid),com_sess_pct,nan(size(suid)),double(regsess),double(waveid)];
                %==================================================1sess=====================23suid======45COM_context1===67COM_context2===89ccfid=====waveid======
            end
        end


        function fh=plot(fc_com_pvsst_stats,reg_com_maps,opt)
            arguments
                fc_com_pvsst_stats
                reg_com_maps
                opt.condense_plot (1,1) logical = false
                opt.omit_reg_wave (1,1) logical = true
                opt.odor_only (1,1) logical = false
            end
            fns={'olf','dur','mixed'};
            idmap=load(fullfile('..','align','reg_ccfid_map.mat'));

            congrusel=pct.su_pairs.get_congru(fc_com_pvsst_stats(:,10:11));
            % incongsel=pct.su_pairs.get_incongru(fc_com_pvsst_stats(:,10:11));
            olf_congru=congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),7:8),2);
            mixed_congru=congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:8),2);
            dur_congru=congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:6),2);
            typesel_mat=[olf_congru,dur_congru,mixed_congru];


            barcnt=[];
            barmm=[];
            barci=[];
            if false %forward and reverse
                for typeIdx=1:3
                    reg_tcom_map=reg_com_maps.(fns{typeIdx});

                    avail_regs=cell2mat(idmap.reg2ccfid.values(reg_tcom_map.keys()));
                    reg_sel=all(ismember(fc_com_pvsst_stats(:,8:9),avail_regs),2);

                    reg_wave_timing=cellfun(@(x) reg_tcom_map(x{1}),...
                        idmap.ccfid2reg.values(...
                        num2cell(fc_com_pvsst_stats(reg_sel,8:9))));
                    fc_com_pvsst_stats(:,6:7)=nan;
                    fc_com_pvsst_stats(reg_sel,6:7)=reg_wave_timing;
                    fini_sel=all(isfinite(fc_com_pvsst_stats(:,4:7)),2);
                    typesel=typesel_mat(:,typeIdx);
                    reg_latency=diff(fc_com_pvsst_stats(fini_sel & typesel,6:7),1,2);
                    fc_latency=diff(fc_com_pvsst_stats(fini_sel & typesel,4:5),1,2);
                    same_reg_sel=fc_com_pvsst_stats(fini_sel & typesel,8)==fc_com_pvsst_stats(fini_sel & typesel,9) & all(fc_com_pvsst_stats(fini_sel & typesel,8:9)>0,2);
                    wave_reg_sel=fc_com_pvsst_stats(fini_sel & typesel,6)<fc_com_pvsst_stats(fini_sel & typesel,7) & all(fc_com_pvsst_stats(fini_sel & typesel,8:9)>0,2);
                    reverse_reg_sel=fc_com_pvsst_stats(fini_sel & typesel,6)>fc_com_pvsst_stats(fini_sel & typesel,7) & all(fc_com_pvsst_stats(fini_sel & typesel,8:9)>0,2);
                    % same region
                    [same_hat,same_ci]=binofit(nnz(fc_latency>0 & same_reg_sel),nnz(same_reg_sel));
                    % reg-wave dir
                    [wave_hat,wave_ci]=binofit(nnz(fc_latency>0 & wave_reg_sel),nnz(wave_reg_sel));
                    % anti-reg-wave dir
                    [rev_hat,rev_ci]=binofit(nnz(fc_latency>0 & reverse_reg_sel),nnz(reverse_reg_sel));
                    % different grouped by regional per-wave timing
                    %         [reg_wave_hat,reg_wave_ci]=binofit(nnz(reg_latency>0 & reverse_reg_sel),nnz(reverse_reg_sel));
                    barcnt=[barcnt;...
                        nnz(fc_latency>0 & same_reg_sel),...
                        nnz(same_reg_sel),...
                        nnz(fc_latency>0 & wave_reg_sel),...
                        nnz(wave_reg_sel)...
                        nnz(fc_latency>0 & reverse_reg_sel),...
                        nnz(reverse_reg_sel)];
                    barmm=[barmm;...
                        same_hat,1-same_hat,...
                        wave_hat,1-wave_hat,...
                        rev_hat,1-rev_hat...
                        ];
                    barci=[barci;same_ci,wave_ci,rev_ci];
                end
            else% same and diff % TODO
                for typeIdx=1:3
                    if opt.odor_only && typeIdx>1
                        continue
                    end
                    reg_tcom_map=reg_com_maps.(fns{typeIdx});
                    avail_regs=cell2mat(idmap.reg2ccfid.values(reg_tcom_map.keys()));
                    reg_sel=all(ismember(fc_com_pvsst_stats(:,8:9),avail_regs),2);
                    reg_wave_timing=cellfun(@(x) reg_tcom_map(x{1}),...
                        idmap.ccfid2reg.values(...
                        num2cell(fc_com_pvsst_stats(reg_sel,8:9))));
                    fc_com_pvsst_stats(:,6:7)=nan; % ensure unique data source
                    fc_com_pvsst_stats(reg_sel,6:7)=reg_wave_timing ;
                    fini_sel=all(isfinite(fc_com_pvsst_stats(:,4:7)),2);
                    typesel=typesel_mat(:,typeIdx);

                    reg_latency=diff(fc_com_pvsst_stats(fini_sel & typesel,6:7),1,2);
                    fc_latency=diff(fc_com_pvsst_stats(fini_sel & typesel,4:5),1,2);
                    % TODO region-averaged latency

                    same_reg_sel=fc_com_pvsst_stats(fini_sel & typesel,8)==fc_com_pvsst_stats(fini_sel & typesel,9) & all(fc_com_pvsst_stats(fini_sel & typesel,8:9)>0,2);
                    diff_reg_sel=fc_com_pvsst_stats(fini_sel & typesel,6)~=fc_com_pvsst_stats(fini_sel & typesel,7) & all(fc_com_pvsst_stats(fini_sel & typesel,8:9)>0,2);

                    % same region
                    [same_hat,same_ci]=binofit(nnz(fc_latency>0 & same_reg_sel),nnz(same_reg_sel));
                    % reg-wave dir
                    [diff_hat,diff_ci]=binofit(nnz(fc_latency>0 & diff_reg_sel),nnz(diff_reg_sel));
                    [diff_reg_hat,diff_reg_ci]=binofit(nnz(reg_latency>0 & diff_reg_sel),nnz(diff_reg_sel));

                    barcnt=[barcnt;...
                        nnz(fc_latency>0 & same_reg_sel),... % 1
                        nnz(same_reg_sel),... % 2
                        nnz(fc_latency>0 & diff_reg_sel),... % 3
                        nnz(diff_reg_sel),... % 4
                        nnz(reg_latency>0 & diff_reg_sel),... % 5
                        nnz(diff_reg_sel)... % 6
                        ];
                    barmm=[barmm;...
                        same_hat,1-same_hat,...
                        diff_hat,1-diff_hat,...
                        diff_reg_hat,1-diff_reg_hat,...
                        ];
                    barci=[barci;same_ci,diff_ci,diff_reg_ci];
                end
            end

            if opt.condense_plot
                fh=figure('Color','w','Position',[100,100,720,235]);
                tiledlayout(1,3);
                nexttile(1,[1,2])
                hold on
                if opt.odor_only
                    bh=bar(1,barmm(1,5:6),'grouped');   % diff region, by reg_tcom
                    errorbar([bh.XEndPoints],[bh.YEndPoints],...
                        [barci(:,5);1-barci(:,5)].'-[bh.YEndPoints],...
                        [barci(:,6);1-barci(:,6)].'-[bh.YEndPoints],'k.')
                    set(gca(),'XTick',1,'XTickLabel',{'Olfactory'},...
                        'YTick',0:0.25:0.75,'YTickLabel',0:25:75)
                else
                    bh=bar(1:3,barmm(:,5:6),'grouped');   % diff region, by reg_tcom
                    errorbar([bh.XEndPoints],[bh.YEndPoints],...
                        [barci(:,5);1-barci(:,5)].'-[bh.YEndPoints],...
                        [barci(:,6);1-barci(:,6)].'-[bh.YEndPoints],'k.')
                    set(gca(),'XTick',1:3,'XTickLabel',{'Olfactory','Duration','Mixed'},...
                        'YTick',0:0.25:0.75,'YTickLabel',0:25:75)
                    xlim([0.5,3.5])
                end
                ylim([0,0.75])
                yline(0.5,'k--')
                
                binocdfp=nan(1,3);
                for typeidx=1:3
                    if typeidx>1 && opt.odor_only
                        continue
                    end
                    binomin=@(x,y) min(barcnt(typeidx,x),barcnt(typeidx,y)-barcnt(typeidx,x));
                    binocdfp(typeidx)=2*binocdf(binomin(5,6),barcnt(typeidx,6),0.5);
                end
                ephys.util.figtable(fh,nexttile(3),binocdfp,'title','binocdf-p')
            else
                fh=figure('Color','w','Position',[100,100,720,235]);
                tiledlayout(1,4)
                %TODO: region-defined panels-> selectivity defined panels
                titles={'Olf.','Dur.','Mixed'};
                for ii=1:3
                    nexttile()
                    hold on
                    bh=bar(1:3,[barmm(ii,[1 3 5]);barmm(ii,[2 4 6])],'grouped');
                    %TODO update
                    errorbar([bh.XEndPoints],[bh.YEndPoints],...
                        [barci(ii,[1 3 5]),1-barci(ii,[2 4 6])]-[bh.YEndPoints],...
                        [barci(ii,[2 4 6]),1-barci(ii,[2 4 6])]-[bh.YEndPoints],'k.')
                    ylim([0,0.75])
                    yline(0.5,'k--')
                    set(gca(),'XTick',1:3,'XTickLabel',{'Within','Cross','Reg-wave'},...
                        'YTick',0:0.25:0.75,'YTickLabel',0:25:75)
                    title(titles{ii});
                    if opt.omit_reg_wave
                        xlim([0.5,2.5])
                    end
                    %         subtitle(sprintf('%d, ',barcnt(:,(1:2)+ii).'));
                end

                disp(num2str(nnz(wave_meta.wave_id>0))+" selective SUs");

                disp("type idx:1, cross reg. condition:1, binocdf each reg.condition:3")
                binocdfp=nan(3,3);
                for typeidx=1:3
                    [~,~,chi2p(typeidx)]=crosstab([zeros(barcnt(typeidx,2),1);...
                        ones(barcnt(typeidx,6),1)],...
                        [(1:barcnt(typeidx,2))>barcnt(typeidx,1),...
                        (1:barcnt(typeidx,4))>barcnt(typeidx,3)]);

                    binomin=@(x,y) min(barcnt(typeidx,x),barcnt(typeidx,y)-barcnt(typeidx,x));
                    binocdfp(typeidx,:)=[2*binocdf(binomin(1,2),barcnt(typeidx,2),0.5),...
                        2*binocdf(binomin(3,4),barcnt(typeidx,4),0.5),...
                        2*binocdf(binomin(5,6),barcnt(typeidx,6),0.5)];
                    %     disp([typeidx,chi2p(typeidx),binocdfp]);
                end
                sgtitle(sprintf('p=%.3f, %.3f, %.3f ',chi2p));
                ephys.util.figtable(fh,nexttile(4),binocdfp,'title','same|fcwave|regwave')
            end
        end

        function alt_plot()
            colors={'r','b','k'};
            titles={'Olfactory','Duration','Multiplexed'};
            % regional wave timing selectivity dependent >>>>>>>>>>>>>>>>
            figure('Color','w','Position',[100,100,1024,235])

            tiledlayout(1,3)
            for typeIdx=1:3
                nexttile();
                hold on

                reg_tcom_map=reg_com_maps.(fns{typeIdx});

                avail_regs=cell2mat(idmap.reg2ccfid.values(reg_tcom_map.keys()));
                reg_sel=all(ismember(fc_com_pvsst_stats(:,8:9),avail_regs),2);

                reg_wave_timing=cellfun(@(x) reg_tcom_map(x{1}),...
                    idmap.ccfid2reg.values(...
                    num2cell(fc_com_pvsst_stats(reg_sel,8:9))));
                fc_com_pvsst_stats(:,6:7)=nan;
                fc_com_pvsst_stats(reg_sel,6:7)=reg_wave_timing;
                fini_sel=all(isfinite(fc_com_pvsst_stats(:,4:7)),2);
                % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                % correlation region wave timing v.s. FC latency

                typesel=typesel_mat(:,typeIdx);

                reg_latency=diff(fc_com_pvsst_stats(fini_sel & typesel,6:7),1,2);
                fc_latency=diff(fc_com_pvsst_stats(fini_sel & typesel,4:5),1,2);
                [N,edges,binidx]=histcounts(reg_latency,-0.35:0.1:0.35);

                fcmm=arrayfun(@(x) mean(fc_latency(binidx==x),'all'),unique(binidx(binidx>0)))./4;
                fcstd=arrayfun(@(x) std(fc_latency(binidx==x)),unique(binidx(binidx>0)))./4;
                fcsem=reshape(fcstd,1,[])./sqrt(N(N>0));
                xep=edges(2:end)-0.05;
                fill([xep(N>0),fliplr(xep(N>0))],[fcmm(:)+fcsem(:);flip(fcmm(:)-fcsem(:))],colors{typeIdx},'EdgeColor','none','FaceAlpha',0.2)
                plot(xep(N>0),fcmm,'Color',colors{typeIdx});
                ylim([-0.8,0.4])
                yline(0,'k--')
                xlabel('Region wave timing latency (s)')
                ylabel('F.C. COM latency (s)')
                set(gca(),'YTick',-0.8:0.4:0.4)
                title(titles{typeIdx});
            end
            sgtitle(num2str(nnz(wave_meta.wave_id>0))+" selective SUs");
        end
        % consistent v. inconsistent
        % congruent v. incongruent
        function alt_plot_2()
            [fwdhat, fwdci]=binofit(barcnt(:,4),sum(barcnt(:,[4 6]),2));
            [revhat, revci]=binofit(barcnt(:,6),sum(barcnt(:,[4 6]),2)); %not really necessary

            bcdfp=[binocdf(min(barcnt(1,4),barcnt(1,6)),sum(barcnt(1,[4,6]),'all'),0.5).*2;...
                binocdf(min(barcnt(2,4),barcnt(2,6)),sum(barcnt(2,[4,6]),'all'),0.5).*2;...
                binocdf(min(barcnt(3,4),barcnt(3,6)),sum(barcnt(3,[4,6]),'all'),0.5).*2];

            %%
            figure() % FC rate along wave direction
            hold on
            % bh=bar([fwdhat,revhat]);
            % errorbar([bh.XEndPoints],[bh.YEndPoints],...
            %     [fwdci(:,1);revci(:,1)].'-[bh.YEndPoints],...
            %     [fwdci(:,2);revci(:,2)].'-[bh.YEndPoints],...
            %     'k.');
            bh=bar(fwdhat,'FaceColor','w','EdgeColor','k');
            errorbar([bh.XEndPoints],[bh.YEndPoints],...
                fwdci(:,1).'-[bh.YEndPoints],...
                fwdci(:,2).'-[bh.YEndPoints],...
                'k.');
            yline(0.5,'k:');
            ylabel('Proportion of F.C (%)')
            set(gca(),'YLim',[0,0.75],'YTick',0:0.25:0.75,'YTickLabel',0:25:75,...
                'XTick',1:3,'XTickLabel',{'Olf.','Dur.','Mixed'})
            % legend(bh,{'Leading- to following reg.','Following- to leading reg.'},...
            %     'Location','northoutside','Orientation','horizontal');
            title(sprintf('%.3f,',bcdfp));
        end
    end
end