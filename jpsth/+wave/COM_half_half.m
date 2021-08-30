% Showcase from s18
if isunix
    if isempty(gcp('nocreate'))
        parpool(50);
    end
    rpts=500;
else
    if isempty(gcp('nocreate'))
        parpool(2)
    end
    rpts=5;
end
to_plot=false;
localshuff=@(x) randsample(x,numel(x));

stats=nan(rpts,2);
[data1hsum,data2hsum]=deal([]);
parfor rpt=1:rpts
    disp(rpt)
    [data1all,data2all]=deal([]);
    com_map=wave.get_com_map('curve',false,'rnd_half',true);
    for sid=1:116
        fs=['s',num2str(sid)];
        s1key=num2cell(intersect(cell2mat(com_map.(fs).s1a.keys),cell2mat(com_map.(fs).s1b.keys)));
        s2key=num2cell(intersect(cell2mat(com_map.(fs).s2a.keys),cell2mat(com_map.(fs).s2b.keys)));
        data1h=[cell2mat(com_map.(fs).s1a.values(s1key)).';cell2mat(com_map.(fs).s2a.values(s2key)).'].*0.25;
        data2h=[cell2mat(com_map.(fs).s1b.values(s1key)).';cell2mat(com_map.(fs).s2b.values(s2key)).'].*0.25;
        data1all=[data1all;data1h];
        data2all=[data2all;data2h];
        if to_plot
            [rd,pd]=corr(data1h,data2h);
            [rs,ps]=corr(data1h,localshuff(data2h));
            fh=figure('Color','w','Position',[32,32,210,210]);
            hold on;
            dh=scatter(data1h,data2h,9,'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
            sh=scatter(data1h,(localshuff(data2h)),9,'o','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
            xlabel('C.O.M. in 1st half trials (s)')
            ylabel('C.O.M. in 2nd half trials (s)')
            legend([dh,sh],{'Real','Shuffle'},'Location','northoutside','Orientation','horizontal')
            text(max(xlim()),max(ylim()),sprintf('%.2f,%.2f,%.2f,%.2f',rd,pd,rs,ps),'HorizontalAlignment','right','VerticalAlignment','top');
            keyboard();
            exportgraphics(fh,'COM_half_half.pdf');
        end
    end
    [rd,pd]=corr(data1all,data2all);
    [rs,ps]=corr(data1all,localshuff(data2all));
    stats(rpt,:)=[rd,rs];
end

