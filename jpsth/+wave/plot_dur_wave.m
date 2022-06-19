function fh=plot_dur_wave(com_map,opt)
arguments
    com_map
    opt.exportgraphics (1,1) logical = false
end
%% global
immat=struct();
immat_anti=struct();
com=struct();
[immat.c1,immat_anti.c1,com.c1,immat.c2,immat_anti.c2,com.c2]=deal([]);
for fsess=reshape(fieldnames(com_map.summed),1,[])
    fs=char(fsess);
    disp(fsess)
    for pref=["c1","c2"]
        if ~isfield(com_map.summed.(fs),pref)
            disp(string(fs)+" missing field "+pref);
            continue;
        end
        curr_key=cell2mat(com_map.summed.(fs).(pref).keys);
        heat=cell2mat(com_map.summed.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
        anti=cell2mat(com_map.summed.(fs).([char(pref),'anticurve']).values(num2cell(curr_key.')));
        COM=cell2mat(values(com_map.summed.(fs).(pref),num2cell(curr_key)));
        %         keyboard()
        immat.(pref)=[immat.(pref);heat];
        immat_anti.(pref)=[immat_anti.(pref);anti];
        com.(pref)=[com.(pref);COM.'];
    end
end

currscale=max(abs([immat.c1(:,1:12),immat_anti.c1(:,1:12);immat.c2(:,1:12),immat_anti.c2(:,1:12)]),[],2);
[immat.c1,immat.c2,immat_anti.c1,immat_anti.c2]...
    =deal(immat.c1./currscale(1:size(immat.c1,1)),immat.c2./currscale(size(immat.c1,1)+1:end),...
    immat_anti.c1./currscale(1:size(immat.c1,1)),immat_anti.c2./currscale(size(immat.c1,1)+1:end));



%% wave
fh=figure('Color','w','Position',[32,32,1400,800]);
[~,sortidx]=sort(com.c1);
plotOne(1,immat.c1(sortidx,:),com.c1(sortidx),'sub_dim',[2,2],'title','FR, preferred','cmap','turbo');
plotOne(2,immat_anti.c1(sortidx,:),com.c1(sortidx),'sub_dim',[2,2],'title','FR, non-preferred','cmap','turbo');
[~,sortidx]=sort(com.c2);
plotOne(4,immat.c2(sortidx,:),com.c2(sortidx),'sub_dim',[2,2],'title','FR, preferred','cmap','turbo');
plotOne(3,immat_anti.c2(sortidx,:),com.c2(sortidx),'sub_dim',[2,2],'title','FR, non-preferred','cmap','turbo');
if opt.exportgraphics
    exportgraphics(gcf(),'dur_wave.pdf','ContentType','vector')
end

end


function plotOne(subidx,imdata,comdata,opt)
arguments
    subidx (1,1) double {mustBeInteger,mustBePositive}
    imdata (:,:) double
    comdata = []
    opt.sub_dim (1,2) = [1,3]
    opt.plot_com (1,1) logical = false
    opt.scale (1,2) double = [-1,1]
    opt.title (1,:) char = []
    opt.cmap (1,:) char = 'turbo'

end
subplot(opt.sub_dim(1),opt.sub_dim(2),subidx);
hold on
colormap(opt.cmap);
gk = fspecial('gaussian', [3 3], 1);
% if size(imdata,2)<20
%     cpos=get(axes,'Position');
%     axes(sh,'Position',[cpos(1),cpos(2),0.5*(diff(cpos([1 3]))),cpos(4)])
% end
imagesc(conv2(imdata,gk,'same'),opt.scale)
if opt.plot_com && exist('comdata','var') && ~isempty(comdata)
    scatter(comdata,1:numel(comdata),2,'o','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
end
if size(imdata,2)>20
    set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
    xline(12.5,'--w','LineWidth',1);
    xlim([0.5,size(imdata,2)+0.5])
else
    xlim([0.5,size(imdata,2)*2+0.5])
    set(gca(),'XTick',[0.5,12.5],'XTickLabel',[0,3]);
end
colorbar();
ylim([0.5,size(imdata,1)+0.5])
if numel(opt.title)>0
    title(opt.title)
end

end