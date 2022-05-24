function out=plot_dur_wave(opt)
arguments
    opt.sess (1,1) double {mustBeMember(opt.sess,1:116)} = 102
    opt.sortby (1,:) char {mustBeMember(opt.sortby,{'3s','6s'})} = '6s'
    opt.exportgraphics (1,1) logical = false
    opt.plot_global (1,1) logical = true
    opt.plot_session (1,1) logical = false
    opt.comb_set (1,:) double {mustBeInteger,mustBePositive} = 1:3
    opt.bootrpt (1,1) double {mustBeInteger,mustBePositive} = 3;
    opt.plot_2d_corr (1,1) logical = true
    opt.plot_corr_dist (1,1) logical = true
end
%% global
out=struct();
out.corr2=struct();
out.corr_dist=struct();

com_map=wave.get_dur_com_map('curve',true);
% tag='Selective in both 3s and 6s trials';
fn='dur_wave.pdf';

immat=struct();
immat_anti=struct();
com=struct();
[immat.d3,immat_anti.d3,com.d3,immat.d6,immat_anti.d6,com.d6]=deal([]);
for fsess=reshape(fieldnames(com_map),1,[])
    fs=char(fsess);
    for pref=["d3","d6"]
        curr_key=cell2mat(com_map.(fs).(pref).keys);
        heat=cell2mat(com_map.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
        anti=cell2mat(com_map.(fs).([char(pref),'anticurve']).values(num2cell(curr_key.')));
        COM=cell2mat(values(com_map.(fs).(pref),num2cell(curr_key)));
        %         keyboard()
        immat.(pref)=[immat.(pref);heat];
        immat_anti.(pref)=[immat_anti.(pref);anti];
        com.(pref)=[com.(pref);COM.'];
    end
end

currscale=max(abs([immat.d3,immat_anti.d3;immat.d6,immat_anti.d6]),[],2);
[immat.d3,immat.d6,immat_anti.d3,immat_anti.d6]...
    =deal(immat.d3./currscale(1:size(immat.d3,1)),immat.d6./currscale(size(immat.d3,1)+1:end),...
    immat_anti.d3./currscale(1:size(immat.d3,1)),immat_anti.d6./currscale(size(immat.d3,1)+1:end));


if opt.plot_global
    %% wave
    fh=figure('Color','w','Position',[32,32,1400,800]);
    [~,sortidx]=sort(com.d3);
    plotOne(1,immat.d3(sortidx,:),com.d3(sortidx),'sub_dim',[2,2],'title','FR, preferred','cmap','turbo');
    plotOne(2,immat_anti.d3(sortidx,:),com.d3(sortidx),'sub_dim',[2,2],'title','FR, non-preferred','cmap','turbo');
    [~,sortidx]=sort(com.d6);
    plotOne(3,immat.d6(sortidx,:),com.d6(sortidx),'sub_dim',[2,2],'title','FR, preferred','cmap','turbo');
    plotOne(4,immat_anti.d6(sortidx,:),com.d6(sortidx),'sub_dim',[2,2],'title','FR, non-preferred','cmap','turbo');    
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