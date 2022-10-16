function fc_com_pvsst_stats=fc_com_reg_wave(wave_meta,com_map,tcom_maps)

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));

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
    for ff=["s1d3","s1d6","s2d3","s2d6","olf_s1","olf_s2","dur_d3","dur_d6"]
        sukeys=com_map.(['s',num2str(sii)]).(ff).com.keys(); % prefered SUid
        susel=ismember(suid,int32(cell2mat(sukeys)));% same dim as suid
        com_sess_pct(susel)=cell2mat(com_map.(['s',num2str(sii)]).(ff).com.values(num2cell(suid(susel)))); % out put is nx2 in dim
    end
    fc_com_pvsst_stats=[fc_com_pvsst_stats;double(sii).*ones(size(suid(:,1))),double(suid),com_sess_pct,nan(size(suid)),double(regsess),double(waveid)];
    %==================================================sess=====================suid=======COM_context1=====COM_context2=====ccfid===========waveid======
end

congrusel=pct.su_pairs.get_congru(fc_com_pvsst_stats(:,10:11));
incongsel=pct.su_pairs.get_incongru(fc_com_pvsst_stats(:,10:11));
mixed_congru=congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:8),2);
olf_congru=congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),7:8),2);
dur_congru=congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:6),2);
typesel_mat=[olf_congru,dur_congru,mixed_congru,incongsel];
tcom_maps_rearr=tcom_maps([2,3,1,2]);

if false
    colors={'r','b','k','m'};
    titles={'Olfactory','Duration','Multiplexed','Incongru-OLF-gradient'};
    % regional wave timing selectivity dependent >>>>>>>>>>>>>>>>
    figure('Color','w','Position',[100,100,1024,235])

    tiledlayout(1,4)
    for typeIdx=1:4
        nexttile();
        hold on

        sel_tcom_map=tcom_maps_rearr{typeIdx};

        avail_regs=cell2mat(idmap.reg2ccfid.values(sel_tcom_map.keys()));
        reg_sel=all(ismember(fc_com_pvsst_stats(:,8:9),avail_regs),2);

        reg_wave_timing=cellfun(@(x) sel_tcom_map(x{1}),...
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
barcnt=[];
barmm=[];
barci=[];
for typeIdx=1:3
    sel_tcom_map=tcom_maps_rearr{typeIdx};

    avail_regs=cell2mat(idmap.reg2ccfid.values(sel_tcom_map.keys()));
    reg_sel=all(ismember(fc_com_pvsst_stats(:,8:9),avail_regs),2);

    reg_wave_timing=cellfun(@(x) sel_tcom_map(x{1}),...
        idmap.ccfid2reg.values(...
        num2cell(fc_com_pvsst_stats(reg_sel,8:9))));
    fc_com_pvsst_stats(:,6:7)=nan;
    fc_com_pvsst_stats(reg_sel,6:7)=reg_wave_timing;
    fini_sel=all(isfinite(fc_com_pvsst_stats(:,4:7)),2);
    typesel=typesel_mat(:,typeIdx);
    %     reg_latency=diff(fc_com_pvsst_stats(fini_sel & typesel,6:7),1,2);
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


figure('Color','w','Position',[100,100,720,235])
tiledlayout(1,3)
%TODO: region-defined panels-> selectivity defined panels
if false
    titles={'Within region','Within region',...
    'Early region to late region','Early region to late region',...
    'Late region to early region','Late region to early region'};
    nexttile()
    hold on
    bh=bar(1:3,barmm(:,(1:2)+ii),'grouped');
    errorbar([bh.XEndPoints],[bh.YEndPoints],...
        [barci(:,1+ii);1-barci(:,1+ii)]-[bh.YEndPoints].',...
        [barci(:,2+ii);1-barci(:,2+ii)]-[bh.YEndPoints].','k.')
    ylim([0,0.8])
    yline(0.5,'k--')
    set(gca(),'XTick',1:3,'XTickLabel',{'OLF','DUR','MIX'},...
        'YTick',0:0.25:0.75,'YTickLabel',0:25:75)
    title(titles{ii+1});
%     subtitle(sprintf('%d, ',barcnt(:,(1:2)+ii).'));
else
    for ii=1:3
        titles={'Olf.','Dur.','Mixed'};
        nexttile()
        hold on
        bh=bar(1:3,[barmm(ii,[1,3,5]);barmm(ii,[2,4,6])],'grouped');
        errorbar([bh.XEndPoints],[bh.YEndPoints],...
            [barci(ii,[1 3 5]),1-barci(ii,[1 3 5])]-[bh.YEndPoints],...
            [barci(ii,[2 4 6]),1-barci(ii,[2 4 6])]-[bh.YEndPoints],'k.')
        ylim([0,0.8])
        yline(0.5,'k--')
        set(gca(),'XTick',1:3,'XTickLabel',{'Within','Wave','Rev.'},...
            'YTick',0:0.25:0.75,'YTickLabel',0:25:75)
        title(titles{ii});
%         subtitle(sprintf('%d, ',barcnt(:,(1:2)+ii).'));
    end
end
disp(num2str(nnz(wave_meta.wave_id>0))+" selective SUs");

disp("type idx:1, cross reg. condition:1, binocdf each reg.condition:3")
for typeidx=1:3
    [~,~,chi2p(typeidx)]=crosstab([zeros(barcnt(typeidx,2),1);...
        ones(barcnt(typeidx,4),1);...
        2.*ones(barcnt(typeidx,6),1)],...
        [(1:barcnt(typeidx,2))>barcnt(typeidx,1),...
        (1:barcnt(typeidx,4))>barcnt(typeidx,3),...
        (1:barcnt(typeidx,6))>barcnt(typeidx,5)]);

    binomin=@(x,y) min(barcnt(typeidx,x),barcnt(typeidx,y)-barcnt(typeidx,x));
    binocdfp=[2*binocdf(binomin(1,2),barcnt(typeidx,2),0.5),...
        2*binocdf(binomin(3,4),barcnt(typeidx,4),0.5),...
        2*binocdf(binomin(5,6),barcnt(typeidx,6),0.5)];
    disp([typeidx,chi2p(typeidx),binocdfp]);
end
sgtitle(sprintf('%.3f',chi2p));

[fwdhat, fwdci]=binofit(barcnt(:,4),sum(barcnt(:,[4 6]),2));
[revhat, revci]=binofit(barcnt(:,6),sum(barcnt(:,[4 6]),2)); %not really necessary

bcdfp=[binocdf(min(barcnt(1,4),barcnt(1,6)),sum(barcnt(1,[4,6]),'all'),0.5).*2;...
    binocdf(min(barcnt(2,4),barcnt(2,6)),sum(barcnt(2,[4,6]),'all'),0.5).*2;...
    binocdf(min(barcnt(3,4),barcnt(3,6)),sum(barcnt(3,[4,6]),'all'),0.5).*2];

%%
figure() % FC rate along wave direction
hold on
bh=bar([fwdhat,revhat]);
errorbar([bh.XEndPoints],[bh.YEndPoints],...
    [fwdci(:,1);revci(:,1)].'-[bh.YEndPoints],...
    [fwdci(:,2);revci(:,2)].'-[bh.YEndPoints],...
    'k.');
yline(0.5,'k:');
ylabel('Proportion of F.C (%)')
set(gca(),'YLim',[0,0.75],'YTick',0:0.25:0.75,'YTickLabel',0:25:75,...
    'XTick',1:3,'XTickLabel',{'Olf.','Dur.','Mixed'})
legend(bh,{'Leading- to following reg.','Following- to leading reg.'},...
    'Location','northoutside','Orientation','horizontal');
title(sprintf('%.3f,',bcdfp));


end

% consistent v. inconsistent

% congruent v. incongruent