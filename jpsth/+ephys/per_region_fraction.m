function collection=per_region_fraction(opt)
arguments
    opt.png (1,1) logical = false
    opt.pdf (1,1) logical = false
    opt.memtype (1,:) char {mustBeMember(opt.memtype,{'any','sust','trans'})} = 'any'
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6 % DPA delay duration
end

persistent collection_ png_ pdf_ memtype_ delay_
if opt.delay==6
    warning('Delay set to default 6')
end

if isempty(collection_) || opt.png~=png_ || opt.pdf ~= pdf_ || ~strcmp(opt.memtype,memtype_) || opt.delay~=delay_
    meta=ephys.util.load_meta('delay',opt.delay);
    fh=figure('Color','w');
    cmap=colormap('jet');
    collection=recurSubregions(meta,1,ismember(meta.reg_tree(1,:),{'BS','CH'}),0,1,cmap,[],opt.memtype);
    blame=vcs.blame();
    save('per_region_fraction_collection.mat','collection','blame');
    colormap('jet')
    ch=colorbar('Ticks',0:0.25:1,'TickLabels',0:0.2:0.8);
    ch.Label.String='Selective fraction';
    sum1=nnz(ismember(meta.reg_tree(1,:),{'CH','BS'}));
    set(gca,'XTick',0.1:0.2:0.9,'XTickLabel',1:5,'YTick',(0:5000:25000)/sum1,'YTickLabel',0:5000:25000)
    xlabel('Allen CCF tree branch level')
    ylabel('Single unit #')
    if ~(opt.png || opt.pdf)
        close(fh)
    end
    if opt.png
        set(fh, 'PaperPosition', [0 0 12 40])    % can be bigger than screen
        print(fh, 'per_region_fraction.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi
    end
    if opt.pdf
        set(gcf, 'PaperSize', [12 40])    % Same, but for PDF output
        print(gcf, 'per_region_fraction.pdf', '-dpdf', '-r300' );   %save file as PDF w/ 300dpi
    end
    collection_=collection;
    memtype_=opt.memtype;
    png_=opt.png;
    pdf_=opt.pdf;
    delay_=opt.delay;
else
    collection=collection_;
end
    function collection=recurSubregions(meta,curr_branch,prev_sel,prev_accu,prev_ratio,cmap,collection,memtype)
        xbound=0:0.2:1;
        if curr_branch==1
            curr_ureg={'BS','CH'};
        else
            curr_ureg=unique(meta.reg_tree(curr_branch,prev_sel));
        end
        switch memtype
            case 'any'
                memset=1:4;
            case 'sust'
                memset=[1,3];
            case 'trans'
                memset=[2,4];
        end
        curr_ratio=cellfun(@(x) nnz(ismember(meta.mem_type(strcmp(meta.reg_tree(curr_branch,:),x)),memset))./...
            nnz(meta.mem_type(strcmp(meta.reg_tree(curr_branch,:),x))>-1),curr_ureg);
            

        [curr_ratio,sidx]=sort(curr_ratio);
        curr_ureg=curr_ureg(sidx);
        curr_count=cellfun(@(x) nnz(meta.mem_type(strcmp(meta.reg_tree(curr_branch,:),x))>-1),curr_ureg);
        collection=[collection;...
            [num2cell(reshape(curr_ratio,[],1)),...
            reshape(curr_ureg,[],1),...
            num2cell(ones(numel(curr_ureg),1)*curr_branch),...
            num2cell(reshape(curr_count,[],1))]];
        curr_accu=prev_accu;
        for curr_idx=1:numel(curr_ureg)
            if isempty(char(curr_ureg(curr_idx))), continue;end
            curr_sel=strcmp(meta.reg_tree(curr_branch,:),curr_ureg(curr_idx));
%             if nnz(curr_sel)<40, continue;end
            sub_ratio=nnz(curr_sel)./nnz(prev_sel)*prev_ratio;
            patch(xbound([curr_branch,curr_branch+1,curr_branch+1,curr_branch]),[curr_accu,curr_accu,curr_accu+sub_ratio,curr_accu+sub_ratio],[-1,-1,-1,-1],cmap(1+floor(curr_ratio(curr_idx)*320),:),'EdgeColor','none');
            text(mean(xbound([curr_branch,curr_branch+1])),curr_accu+sub_ratio/2,1,curr_ureg(curr_idx),'HorizontalAlignment','left','VerticalAlignment','middle')
            if curr_branch<5
                collection=recurSubregions(meta,curr_branch+1,curr_sel,curr_accu,sub_ratio,cmap,collection,memtype);
            end
            curr_accu=curr_accu+sub_ratio;
        end
    end
end
