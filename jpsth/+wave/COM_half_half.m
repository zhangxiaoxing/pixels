function COM_half_half(opt)
arguments
    opt.new_data (1,1) logical =true
    opt.plot (1,1) logical = true
    opt.corr_curve (1,1) logical = false
end
if opt.new_data
    stats=gendata();
    if isunix
        save('COM_half_half.mat','stats');
    end
end


if opt.plot %isfile('COM_half_half.mat')
    load('COM_half_half.mat','stats');
    %TODO else generate small local dataset for showcase?
    sessmm=cell2mat(arrayfun(@(x) mean(stats(stats(:,1)==x,3:8)),unique(stats(:,1)),'UniformOutput',false));
%     mm=mean(sessmm);
%     sem=std(sessmm);
    fh=figure('Color','w','Position',[32,32,420,100]);
    for idx=1:3:4
        subplot(1,2,(idx>1)+1);
        hold on
%         bar(mm(idx:idx+2),0.7,'FaceColor','w','EdgeColor','k','LineWidth',1);
%         errorbar(1:3,mm(idx:idx+2),sem(idx:idx+2),'k.','CapSize',20,'Color',[0.5,0.5,0.5])
        boxplot(sessmm(:,idx:idx+2));
        if idx>1
            ylabel('FR Pearson''s r')
        else
            ylabel('TCOM Pearson''s r')
        end
        set(gca(),'XTick',1:3,'XTickLabel',{'Correct','Error','Shuffle'});
        yline(0,'-k')
        xlim([0.5,3.5]);
        ylim([-0.1,1])
    end
    exportgraphics(fh,'COM_half_half_stats.pdf')
end
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
stats=[];
% [data1hsum,data2hsum]=deal([]);
for rpt=1:rpts
    disp(rpt)
    [data1hall,data2hall,dataeall]=deal([]);
    com_map=wave.get_com_map('curve',true,'rnd_half',true,'wave','any6','delay',6);
    for fn=reshape(fieldnames(com_map),1,[])
        fs=fn{1};
        s1key=num2cell(intersect(cell2mat(com_map.(fs).s1a.keys),intersect(cell2mat(com_map.(fs).s1b.keys),cell2mat(com_map.(fs).s1e.keys))));
        s2key=num2cell(intersect(cell2mat(com_map.(fs).s2a.keys),intersect(cell2mat(com_map.(fs).s2b.keys),cell2mat(com_map.(fs).s2e.keys))));
        if isempty(s1key) || isempty(s2key)
            continue
        end
        data1h=[cell2mat(com_map.(fs).s1a.values(s1key)).';cell2mat(com_map.(fs).s2a.values(s2key)).'].*0.25;
        data2h=[cell2mat(com_map.(fs).s1b.values(s1key)).';cell2mat(com_map.(fs).s2b.values(s2key)).'].*0.25;
        datae=[cell2mat(com_map.(fs).s1e.values(s1key)).';cell2mat(com_map.(fs).s2e.values(s2key)).'].*0.25;

        data1hall=[data1hall;data1h];
        data2hall=[data2hall;data2h];
        dataeall=[dataeall;datae];


        cdata1h=[cell2mat(com_map.(fs).s1acurve.values(s1key.'));cell2mat(com_map.(fs).s2acurve.values(s2key.'))];
        cdata2h=[cell2mat(com_map.(fs).s1bcurve.values(s1key.'));cell2mat(com_map.(fs).s2bcurve.values(s2key.'))];
        cdatae=[cell2mat(com_map.(fs).s1ecurve.values(s1key.'));cell2mat(com_map.(fs).s2ecurve.values(s2key.'))];
        if isempty(cdata1h) || isempty(cdata2h) || isempty(cdatae)
            keyboard()
        end

        scaler=max(abs([cdata1h,cdata2h,cdatae]),[],2);

        % curve corr data
%         keyboard()
        cdata1h=cdata1h./scaler;
        cdata2h=cdata2h./scaler;
        cdatae=cdatae./scaler;

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
        % TODO stats over session ?
        rd=corr(data1hall,data2hall);
        rs=corr(data1hall,localshuff(data2hall));
        re=corr(data1hall,dataeall);

        rcd=corr(cdata1h(:),cdata2h(:));
        cdshuf=cdata2h(localshuff(1:size(cdata2h,1)),:);
        rcs=corr(cdata1h(:),cdshuf(:));
        rce=corr(cdata1h(:),cdatae(:));
        
        sess=str2double(replace(fs,'s',''));
        
        stats=[stats;sess, rpt,rd,re,rs,rcd,rce,rcs];
        assignin('base','stats',stats);

    end

    
end
end