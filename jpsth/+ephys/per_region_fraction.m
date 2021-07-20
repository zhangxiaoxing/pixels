% meta=ephys.util.load_meta();
% ureg1=["BS","CH"];
% ratio1=arrayfun(@(x) nnz(meta.mem_type(strcmp(meta.reg_tree(1,:),x))>0)./...
%     nnz(meta.mem_type(strcmp(meta.reg_tree(1,:),x))>-1),ureg1);
% [ratio1,sidx]=sort(ratio1);
% ureg1=ureg1(sidx);
%
% sum1=sum(ismember(meta.reg_tree(1,:),ureg1));
% xbound=0:0.2:1;

cmap=colormap('jet');
fh=figure('Color','w');
accu1=0;
for r1idx=1:numel(ureg1)
    
    r1sel=strcmp(meta.reg_tree(1,:),ureg1(r1idx));
    cratio1=nnz(r1sel)./sum1;
    patch(xbound([1,2,2,1]),[accu1,accu1,accu1+cratio1,accu1+cratio1],cmap(1+floor(ratio1(r1idx)*320),:),'EdgeColor','none');
    text(mean(xbound([1,2])),accu1+cratio1/2,ureg1(r1idx),'HorizontalAlignment','center','VerticalAlignment','middle')
    %TODO tier2
    ureg2=unique(meta.reg_tree(2,r1sel));
    ratio2=arrayfun(@(x) nnz(meta.mem_type(strcmp(meta.reg_tree(2,:),x))>0)./...
        nnz(meta.mem_type(strcmp(meta.reg_tree(2,:),x))>-1),ureg2);
    [ratio2,sidx]=sort(ratio2);
    ureg2=ureg2(sidx);
    accu2=accu1;
    for r2idx=1:numel(ureg2)
        r2sel=strcmp(meta.reg_tree(2,:),ureg2(r2idx));
        cratio2=nnz(r2sel)./nnz(r1sel)*cratio1;
        patch(xbound([2,3,3,2]),[accu2,accu2,accu2+cratio2,accu2+cratio2],cmap(1+floor(ratio2(r2idx)*320),:),'EdgeColor','none');
        text(mean(xbound([2,3])),accu2+cratio2/2,ureg2(r2idx),'HorizontalAlignment','center','VerticalAlignment','middle')
        %TODO tier3
        ureg3=unique(meta.reg_tree(3,r2sel));
        ratio3=arrayfun(@(x) nnz(meta.mem_type(strcmp(meta.reg_tree(3,:),x))>0)./...
            nnz(meta.mem_type(strcmp(meta.reg_tree(3,:),x))>-1),ureg3);
        [ratio3,sidx]=sort(ratio3);
        ureg3=ureg3(sidx);
        accu3=accu2;
        for r3idx=1:numel(ureg3)
            if isempty(char(ureg3(r3idx))), continue;end
            r3sel=strcmp(meta.reg_tree(3,:),ureg3(r3idx));
            cratio3=nnz(r3sel)./nnz(r2sel)*cratio2;
            patch(xbound([3,4,4,3]),[accu3,accu3,accu3+cratio3,accu3+cratio3],cmap(1+floor(ratio3(r3idx)*320),:),'EdgeColor','none');
            text(mean(xbound([3,4])),accu3+cratio3/2,ureg3(r3idx),'HorizontalAlignment','center','VerticalAlignment','middle')
            %TODO tier4
            ureg4=unique(meta.reg_tree(4,r3sel));
            ratio4=arrayfun(@(x) nnz(meta.mem_type(strcmp(meta.reg_tree(4,:),x))>0)./...
                nnz(meta.mem_type(strcmp(meta.reg_tree(4,:),x))>-1),ureg4);
            [ratio4,sidx]=sort(ratio4);
            ureg4=ureg4(sidx);
            accu4=accu3;
            for r4idx=1:numel(ureg4)
                if isempty(char(ureg4(r4idx))), continue;end
                r4sel=strcmp(meta.reg_tree(4,:),ureg4(r4idx));
                cratio4=nnz(r4sel)./nnz(r3sel)*cratio3;
                patch(xbound([4,5,5,4]),[accu4,accu4,accu4+cratio4,accu4+cratio4],[-1,-1,-1,-1],cmap(1+floor(ratio4(r4idx)*320),:),'EdgeColor','none');
                text(min(xbound([4,5]))+rem(r4idx,4)*0.05,accu4+cratio4/2,ureg4(r4idx),'HorizontalAlignment','left','VerticalAlignment','middle')
                %TODO tier5
                ureg5=unique(meta.reg_tree(5,r4sel));
                ratio5=arrayfun(@(x) nnz(meta.mem_type(strcmp(meta.reg_tree(5,:),x))>0)./...
                    nnz(meta.mem_type(strcmp(meta.reg_tree(5,:),x))>-1),ureg5);
                [ratio5,sidx]=sort(ratio5);
                ureg5=ureg5(sidx);
                accu5=accu4;
                for r5idx=1:numel(ureg5)
                    if isempty(char(ureg5(r5idx))), continue;end
                    r5sel=strcmp(meta.reg_tree(5,:),ureg5(r5idx));
                    cratio5=nnz(r5sel)./nnz(r4sel)*cratio4;
                    patch(xbound([5,6,6,5]),[accu5,accu5,accu5+cratio5,accu5+cratio5],[-1,-1,-1,-1],cmap(1+floor(ratio5(r5idx)*320),:),'EdgeColor','none');
                    text(min(xbound([5,6]))+rem(r5idx,4)*0.05,accu5+cratio5/2,ureg5(r5idx),'HorizontalAlignment','left','VerticalAlignment','middle')
                    accu5=accu5+cratio5;
                end
                accu4=accu4+cratio4;
            end
            accu3=accu3+cratio3;
        end
        accu2=accu2+cratio2;
    end
    accu1=accu1+cratio1;
end
colormap('jet')
ch=colorbar('Ticks',0:0.25:1,'TickLabels',0:0.2:0.8);
ch.Label.String='Selective fraction';
set(gca,'XTick',0.1:0.2:0.9,'XTickLabel',1:5,'YTick',(0:5000:25000)/sum1,'YTickLabel',0:5000:25000)
xlabel('Allen CCF tree branch level')
ylabel('Single unit #')


set(fh, 'PaperPosition', [0 0 12 40])    % can be bigger than screen 
% set(gcf, 'PaperSize', [10 30])    % Same, but for PDF output
print(fh, 'fraction.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi
% print(gcf, 'fraction.pdf', '-dpdf', '-r300' );   %save file as PDF w/ 300dpi
