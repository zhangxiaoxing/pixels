classdef fc_com_reg_wave < handle
    methods (Static)
        function fc_reg_tcom=stats(wrs_mux_meta,reg_com_maps,opt)
            arguments
                wrs_mux_meta
                reg_com_maps
                opt.delay (1,1) double {mustBeMember(opt.delay,[3 6])} = 6
                opt.odor_only (1,1) logical = false
                opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT'})} = 'WT'
            end

            idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
            [sig,~]=bz.load_sig_sums_conn_file("criteria",opt.criteria);
            sig=bz.join_fc_waveid(sig,wrs_mux_meta.wave_id,"criteria",opt.criteria);
                
            grey=reg_com_maps.("tcom"+opt.delay+"_maps").odor_only.keys;

            tokeep=all(sig.waveid==0,2) ...
                | all(ismember(sig.waveid,[1 5]),2)...
                | all(ismember(sig.waveid,[2 5]),2)...
                | all(ismember(sig.waveid,[3 6]),2)...
                | all(ismember(sig.waveid,[4 6]),2);
            
            assert(opt.odor_only,"unfinished")

            fc_reg_tcom=[];
            for fcidx=reshape(find(tokeep),1,[])
                if any(sig.waveid(fcidx,:)==0) && any(sig.waveid(fcidx,:)>0)
                    continue
                end
                if any(sig.reg(fcidx,5,:)==0,'all')
                    continue;
                end
                regs=cellfun(@(x) x{1},idmap.ccfid2reg.values(num2cell(squeeze(sig.reg(fcidx,5,:)))),'UniformOutput',false);
                if ~all(ismember(regs,grey),'all')
                    continue
                end
                reg_tcom=cell2mat(reg_com_maps.("tcom"+opt.delay+"_maps").odor_only.values(regs));
                row=cell2table({sig.sess(fcidx),sig.suid(fcidx,:),sig.waveid(fcidx,:),regs.',reg_tcom.'},...
                    'VariableNames',{'sess','cid','waveid','reg','reg_tcom'});
                fc_reg_tcom=[fc_reg_tcom;row];
            end
        end

        function [barmm,barci,barcnt]=sums(fc_reg_tcom,opt)
            arguments
                fc_reg_tcom
                opt.odor_only (1,1) logical = false
            end
            congru_sel=all(fc_reg_tcom.waveid>0,2);
            nm_sel=~congru_sel;

            cross_sel=~strcmp(fc_reg_tcom.reg(:,1),fc_reg_tcom.reg(:,2));

            consist_sel=diff(fc_reg_tcom.reg_tcom,1,2)>0;
            inconsist_sel=diff(fc_reg_tcom.reg_tcom,1,2)<0;

            assert(opt.odor_only,"unfinished");
            
            barcnt=[];
            barmm=[];
            barci=[];

            congru_consistent=congru_sel & cross_sel & consist_sel;
            congru_inconsistent=congru_sel & cross_sel & inconsist_sel;

            nonmem_consistent=nm_sel & cross_sel & consist_sel;
            nonmem_inconsistent=nm_sel & cross_sel & inconsist_sel;

            [chat,cci]=binofit(nnz(congru_consistent),nnz(congru_consistent | congru_inconsistent));
            [cihat,cici]=binofit(nnz(congru_inconsistent),nnz(congru_consistent | congru_inconsistent));

            [nhat,nci]=binofit(nnz(nonmem_consistent),nnz(nonmem_consistent | nonmem_inconsistent));
            [nihat,nici]=binofit(nnz(nonmem_inconsistent),nnz(nonmem_consistent | nonmem_inconsistent));

            barmm=[chat,cihat;nhat,nihat];
            barci=[cci;cici;nci;nici];
            barcnt=[nnz(congru_consistent),nnz(congru_consistent | congru_inconsistent);...
                nnz(congru_inconsistent),nnz(congru_consistent | congru_inconsistent);...
                nnz(nonmem_consistent),nnz(nonmem_consistent | nonmem_inconsistent);...
                nnz(nonmem_inconsistent),nnz(nonmem_consistent | nonmem_inconsistent)];

        end
        
        function fh=plot(barmm,barci,barcnt)
            arguments
                barmm
                barci
                barcnt
            end


            fh=figure('Color','w','Position',[100,100,720,235]);
            tiledlayout(1,3);
            nexttile(1,[1,2])
            hold on

            bh=bar(barmm,'grouped');   % diff region, by reg_tcom
            errorbar([bh.XEndPoints],[bh.YEndPoints],...
                barci([1 3 2 4],1).'-[bh.YEndPoints],...
                barci([1 3 2 4],2).'-[bh.YEndPoints],'k.')
            set(gca(),'XTick',1:2,'XTickLabel',{'Memory','Nonmemory'},...
                'YTick',0:0.25:0.75,'YTickLabel',0:25:75,'YLim',[0,0.75])

            yline(0.5,'k--')

            [~,~,pcongru]=crosstab(reshape(repmat([0,1],barcnt(1,2),1),[],1),[1:barcnt(1,2)>barcnt(1,1),1:barcnt(2,2)>barcnt(2,1)]);
            [~,~,pnonmem]=crosstab(reshape(repmat([0,1],barcnt(3,2),1),[],1),[1:barcnt(3,2)>barcnt(3,1),1:barcnt(4,2)>barcnt(4,1)]);

            [~,~,pgroup]=crosstab([zeros(barcnt(1,2),1);ones(barcnt(3,2),1)],[1:barcnt(1,2)>barcnt(1,1),1:barcnt(3,2)>barcnt(3,1)]);

            title(sprintf('con%.4f,nonmem%.4f,group%.4f',pcongru,pnonmem,pgroup))


        end

        function alt_plot()
            colors={'r','b','k'};
            titles={'Olfactory','Duration','Multiplexed'};
            % regional wave timing selectivity dependent >>>>>>>>>>>>>>>>
            figure('Color','w','Position',[100,100,1024,235])

            tiledlayout(1,3)
            for typeIdx=1:3
                nexttile()
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


        function fh=run_all(opt)
            arguments
                opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT'})} = 'WT'
            end
            switch opt.criteria
                case 'WT'
                    load(fullfile('binary','wrs_mux_meta.mat'),'wrs_mux_meta');
                    sel_meta=wrs_mux_meta;
                    sfn=fullfile('binary','SC_consistent_inconsistent.fig');
                case 'Learning'
                    sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);
                    sfn=fullfile('binary','LN_SC_consistent_inconsistent.fig');
                otherwise
                    keyboard();
            end

            reg_com_maps=wave.get_reg_com_maps(sel_meta,"criteria",opt.criteria);
            scstats6=fc.fc_com_reg_wave.stats(sel_meta,reg_com_maps,'delay',6,'odor_only',true,'criteria',opt.criteria);
            scstats3=fc.fc_com_reg_wave.stats(sel_meta,reg_com_maps,'delay',3,'odor_only',true,'criteria',opt.criteria);
            scstats=[scstats6;scstats3];

            [barmm,barci,barcnt]=fc.fc_com_reg_wave.sums(scstats,"odor_only",true);
            fh=fc.fc_com_reg_wave.plot(barmm,barci,barcnt);
            
            savefig(fh,sfn);
            if false
                types=repmat("Undefined",size(scstats,1),1);
                types(pct.su_pairs.get_congru(scstats.waveid))="Congruent";
                types(pct.su_pairs.get_incongru(scstats.waveid))="Incongruent";
                types(pct.su_pairs.get_nonmem(scstats.waveid))="Nonmemory";

                scstats=[scstats(:,1:2),array2table(types,'VariableNames',{'Type'}),scstats(:,4:5)];

                fid=fopen(fullfile('binary','upload','F2G_SC_consistent_inconsistent.json'),'w');
                fprintf(fid,jsonencode(scstats));
                fclose(fid);

            end

        end


    end
end