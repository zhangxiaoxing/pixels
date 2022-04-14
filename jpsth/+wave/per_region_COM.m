% TODO: ???
function [collection,com_meta]=per_region_COM(opt)
arguments
    opt.keep_figure (1,1) logical = false
    opt.pdf (1,1) logical = false
    opt.png (1,1) logical = true
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.stats_method (1,:) char {mustBeMember(opt.stats_method,{'mean','median'})} = 'mean';
    opt.selidx (1,1) logical = false % calculate COM of selectivity index
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6 % COM in trials of delay
    opt.wave (1,:) char {mustBeMember(opt.wave,{'both','only3','only6','any3','any6','any'})}  % su of wave
end

persistent com_meta_ collection_ opt_
if opt.delay==6
    warning('Delay set to default 6')
end


if isempty(com_meta_) || isempty(collection_) || ~isequaln(opt,opt_)
    
    % per region COM
    com_map=wave.get_com_map('per_sec_stats',false,'decision',opt.decision,'selidx',opt.selidx,'delay',opt.delay,'wave',opt.wave);
    % COM->SU->region
    meta=ephys.util.load_meta();
    meta.sessid=cellfun(@(x) ephys.path2sessid(x),meta.allpath);
    sess=fieldnames(com_map);
    
    com_meta=cell(0,9);
    for si=1:numel(sess)
        sid=str2double(replace(sess{si},'s',''));
        allcid=meta.allcid(meta.sessid==sid);
        allreg=meta.reg_tree(:,meta.sessid==sid);

        s1id=cell2mat(com_map.(sess{si}).s1.keys()).';
        s1com=cell2mat(com_map.(sess{si}).s1.values()).';
        s2id=cell2mat(com_map.(sess{si}).s2.keys()).';
        s2com=cell2mat(com_map.(sess{si}).s2.values()).';
        
        [~,locs1]=ismember(uint16(s1id),allcid);
        [~,locs2]=ismember(uint16(s2id),allcid);

        com_meta=[com_meta;num2cell([sid*ones(size(s1id)),double(s1id),s1com]),allreg(:,locs1).'];
        com_meta=[com_meta;num2cell([sid*ones(size(s2id)),double(s2id),s2com]),allreg(:,locs2).'];
    end
    
    
    fh=figure('Color','w');
    cmap=colormap('turbo');
    CHratio=nnz(strcmp(com_meta(:,4),'CH'))./size(com_meta,1);
    collectionCH=recurSubregions('CH',strcmp(com_meta(:,4),'CH'),0,CHratio,cmap,com_meta,[],opt);
    collection=recurSubregions('BS',strcmp(com_meta(:,4),'BS'),CHratio,1-CHratio,cmap,com_meta,collectionCH,opt);
    
    blame=vcs.blame();
    save('per_region_com_collection.mat','collection','blame');
    colormap('turbo')
    ch=colorbar('Ticks',0:0.25:1,'TickLabels',1:5);
    ch.Label.String='Mean F.R. COM';
    sum1=nnz(ismember(com_meta(:,4),{'CH','BS'}));
    set(gca,'XTick',0.1:0.2:0.9,'XTickLabel',1:5,'YTick',(0:2000:8000)/sum1,'YTickLabel',0:2000:8000,'YDir','reverse')
    xlabel('Allen CCF tree branch level')
    ylabel('Single unit #')
    
    if opt.png
        set(fh, 'PaperPosition', [0 0 12 40])    % can be bigger than screen
        print(fh, 'per_region_COM.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi
    elseif opt.pdf
        set(gcf, 'PaperSize', [12 40])    % Same, but for PDF output
        print(gcf, 'per_region_COM.pdf', '-dpdf', '-r300' );   %save file as PDF w/ 300dpi
    end
    if ~opt.keep_figure
        close(fh)
    end
    collection_=collection;
    com_meta_=com_meta;
    opt_=opt;
else
    collection=collection_;
    com_meta=com_meta_;
end
end


function collection=recurSubregions(curr_branch,prev_sel,prev_accu,prev_ratio,cmap,com_meta,collection,opt)
% arguments
%     curr_branch
%     prev_sel
%     prev_accu
%     prev_ratio
%     cmap
%     com_meta
%     collection
%     opt.stats_method (1,:) char {mustBeMember(opt.stats_method,{'mean','median'})} = 'mean';
% end
    
xbound=0:0.2:2;
if ~isnumeric(curr_branch) && ismember(curr_branch,{'CH','BS'})
    curr_ureg={curr_branch};
    curr_branch=1;
else
    curr_ureg=unique(com_meta(prev_sel,curr_branch+3));
end
if strcmp(opt.stats_method,'mean')
    curr_com_all=cellfun(@(x) mean([com_meta{strcmp(com_meta(:,curr_branch+3),x),3}]),curr_ureg);
else
    curr_com_all=cellfun(@(x) median([com_meta{strcmp(com_meta(:,curr_branch+3),x),3}]),curr_ureg);    
end
[curr_com_all,sidx]=sort(curr_com_all);
curr_ureg=curr_ureg(sidx);
curr_count=cellfun(@(x) nnz(strcmp(com_meta(:,curr_branch+3),x)),curr_ureg);
collection=[collection;...
    [num2cell(reshape(curr_com_all,[],1)),...
    reshape(curr_ureg,[],1),...
    num2cell(ones(numel(curr_ureg),1)*curr_branch),...
    num2cell(reshape(curr_count,[],1))]];
curr_accu=prev_accu;
for curr_idx=1:numel(curr_ureg)
    if isempty(char(curr_ureg(curr_idx))), continue;end
    curr_sel=strcmp(com_meta(:,curr_branch+3),curr_ureg(curr_idx));
%     if nnz(curr_sel)<20, continue;end
    sub_ratio=nnz(curr_sel)./nnz(prev_sel)*prev_ratio;
    patch(xbound([curr_branch,curr_branch+1,curr_branch+1,curr_branch]),[curr_accu,curr_accu,curr_accu+sub_ratio,curr_accu+sub_ratio],[-1,-1,-1,-1],cmap(com2cmapidx(curr_com_all(curr_idx),opt),:),'EdgeColor','none');
    %     if curr_branch>3
    %         text(min(xbound([curr_branch,curr_branch+1]))+rem(curr_idx,4)*0.05,curr_accu+sub_ratio/2,1,curr_ureg(curr_idx),'HorizontalAlignment','left','VerticalAlignment','middle')
    %     else
    text(mean(xbound([curr_branch,curr_branch+1])),curr_accu+sub_ratio/2,1,curr_ureg(curr_idx),'HorizontalAlignment','left','VerticalAlignment','middle')
    %     end
    if curr_branch<5
        collection=recurSubregions(curr_branch+1,curr_sel,curr_accu,sub_ratio,cmap,com_meta,collection,opt);
    end
    curr_accu=curr_accu+sub_ratio;
end
end

function out=com2cmapidx(com,opt)
if opt.decision
    out=floor(com/4*255)+1;
else
    out=floor((com-4)/16*255)+1;
end
if out<1,out=1;end
if out>256,out=256;end
end