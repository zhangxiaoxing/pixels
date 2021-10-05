if new_data
    stats=gendata();
    save('COM_half_half.mat','stats');
end

% Showcase from s18
if isfile('COM_half_half.mat')
    load('COM_half_half.mat','stats');
    %TODO else generate small local dataset for showcase?
    mm=mean(stats);
    ci=bootci(500,@(x) mean(x), stats);
    fh=figure('Color','w','Position',[32,32,210,60]);
    hold on
    bar(mm,0.7,'FaceColor','w','EdgeColor','k','LineWidth',1);
    errorbar(1:3,mm,ci(1,:)-mm,ci(2,:)-mm,'ko','CapSize',20)
    ylabel('Pearson r')
    set(gca(),'XTick',1:3,'XTickLabel',{'Correct','Error','Shuffle'});
    xlim([0.5,3.5]);
    exportgraphics(fh,'COM_half_half_stats.pdf')
end

function stats=gendata(opt)
arguments
    opt.to_plot (1,1) logical = false
end

if isunix
    if isempty(gcp('nocreate'))
        parpool(50);
    end
    rpts=500;
else
    %     if isempty(gcp('nocreate'))
    %         parpool(2)
    %     end
    rpts=5;
end

localshuff=@(x) randsample(x,numel(x));
stats=nan(rpts,3);
% [data1hsum,data2hsum]=deal([]);
parfor rpt=1:rpts
    disp(rpt)
    [data1hall,data2hall,dataeall]=deal([]);
    com_map=wave.get_com_map('curve',false,'rnd_half',true);
    for sid=1:116
        fs=['s',num2str(sid)];
        s1key=num2cell(intersect(cell2mat(com_map.(fs).s1a.keys),intersect(cell2mat(com_map.(fs).s1b.keys),cell2mat(com_map.(fs).s1e.keys))));
        s2key=num2cell(intersect(cell2mat(com_map.(fs).s2a.keys),intersect(cell2mat(com_map.(fs).s2b.keys),cell2mat(com_map.(fs).s2e.keys))));
        data1h=[cell2mat(com_map.(fs).s1a.values(s1key)).';cell2mat(com_map.(fs).s2a.values(s2key)).'].*0.25;
        data2h=[cell2mat(com_map.(fs).s1b.values(s1key)).';cell2mat(com_map.(fs).s2b.values(s2key)).'].*0.25;
        datae=[cell2mat(com_map.(fs).s1e.values(s1key)).';cell2mat(com_map.(fs).s2e.values(s2key)).'].*0.25;
        
        data1hall=[data1hall;data1h];
        data2hall=[data2hall;data2h];
        dataeall=[dataeall;datae];
        if opt.to_plot
            [rd,pd]=corr(data1h,data2h);
            [rs,ps]=corr(data1h,localshuff(data2h));
            [re,pe]=corr(data1h,datae);
            fh=figure('Color','w','Position',[32,32,210,210]);
            hold on;
            dh=scatter(data1h,data2h,9,'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
            sh=scatter(data1h,(localshuff(data2h)),9,'o','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
            eh=scatter(data1h,datae,9,'o','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
            xlabel('C.O.M. in 1st half trials (s)')
            ylabel('C.O.M. in 2nd half trials (s)')
            legend([dh,eh,sh],{'Correct','Error','Shuffle'},'Location','northoutside','Orientation','horizontal')
            text(max(xlim()),max(ylim()),sprintf('%.2f,%.2f,%.2f,%.2f,%.2f,%.2f',rd,pd,re,pe,rs,ps),'HorizontalAlignment','right','VerticalAlignment','top');
            keyboard();
            exportgraphics(fh,'COM_half_half.pdf');
        end
    end
    [rd,pd]=corr(data1hall,data2hall);
    [rs,ps]=corr(data1hall,localshuff(data2hall));
    [re,pe]=corr(data1hall,dataeall);
    stats(rpt,:)=[rd,re,rs];
end
end