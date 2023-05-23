
function fh=plot_pct_wave(pct_meta,com_map,opt)
arguments
%     eff_meta
    pct_meta
    com_map
    opt.comb_set (1,:) double {mustBeInteger,mustBePositive} = 1
    opt.flex_sort (1,1) logical = true
    opt.xlim (1,1) double {mustBeMember(opt.xlim,[3 6])} = 6
    opt.lesser_grade (1,1) logical = false
    opt.scale (1,2) double = [-1,1]
    opt.gauss2d (1,1) logical = false

end
%% global
fss=reshape(fieldnames(com_map),1,[]);
for plot_id=opt.comb_set
    imdata=struct();
    switch plot_id
        case 1 % multiplexed
            fh.("wave"+string(plot_id))=figure('Color','w','Position',[32,32,1200,400]);
            tiledlayout(4,4)
            tags=["s1d3","s1d6","s2d3","s2d6"];
            for prefidx=1:4
                pref=tags(prefidx);
                imdata.(pref)=struct();
                imdata.(pref).com=[];
                for ff=["s1d3","s1d6","s2d3","s2d6"]
                    imdata.(pref).(ff)=[];
                end
                for fs1=fss
                    fs=char(fs1);
                    if ~isfield(com_map.(fs),pref)
                        continue
                    end
                    if opt.flex_sort 
                        curr_key=com_map.(fs).(pref).com4plot.keys();
                        imdata.(pref).com=cat(1,imdata.(pref).com,cell2mat(values(com_map.(fs).(pref).com4plot,curr_key.')));
                    else
                        curr_key=com_map.(fs).(pref).com.keys();
                        imdata.(pref).com=cat(1,imdata.(pref).com,cell2mat(values(com_map.(fs).(pref).com,curr_key.')));
                    end
                    for ff=["s1d3","s1d6","s2d3","s2d6"]
                        imdata.(pref).(ff)=cat(1,imdata.(pref).(ff),cell2mat(com_map.(fs).(pref).(ff).values(curr_key.')));
                    end
                end

                [~,sortidx]=sort(imdata.(pref).com);
                for ffidx=1:4
                    nexttile((prefidx-1)*4+ffidx);
                    if opt.xlim==3
                        plotOne(imdata.(pref).(tags(ffidx))(sortidx,1:12),'xlim',3,'scale',opt.scale,'gauss2d',opt.gauss2d);
                    else
                        plotOne(imdata.(pref).(tags(ffidx))(sortidx,:),'scale',opt.scale,'gauss2d',opt.gauss2d);
                    end
                end
            end
        case 2 % olf
            fh.("wave"+string(plot_id))=figure('Color','w','Position',[32,32,1200,400]);
            tiledlayout(4,4)
            tags=["olf_s1","olf_s2","dur_d3","dur_d6"];
            for prefidx=1:4
                pref=tags(prefidx);
                imdata.(pref)=struct();
                imdata.(pref).com=[];
                imdata.(pref).com6=[];
                for ff=["s1d3","s1d6","s2d3","s2d6"]
                    imdata.(pref).(ff)=[];
                end
                for fs1=fss
                    fs=char(fs1);
                    if ~isfield(com_map.(fs),pref)
                        continue
                    end
                    if opt.flex_sort && startsWith(pref,"olf")
                        curr_key=num2cell(intersect(cell2mat(com_map.(fs).(pref).com3.keys()),...
                            cell2mat(com_map.(fs).(pref).com6.keys())));
                        imdata.(pref).com=cat(1,imdata.(pref).com,...
                            cell2mat(values(com_map.(fs).(pref).com3,curr_key.')));
%                         keyboard()
                        imdata.(pref).com6=cat(1,imdata.(pref).com6,...
                            cell2mat(values(com_map.(fs).(pref).com6,curr_key.')));
                    elseif opt.flex_sort
                        curr_key=com_map.(fs).(pref).com.keys();
                        imdata.(pref).com=cat(1,imdata.(pref).com,...
                            cell2mat(values(com_map.(fs).(pref).com4plot,curr_key.')));
                    else
                        curr_key=com_map.(fs).(pref).com.keys();
                        imdata.(pref).com=cat(1,imdata.(pref).com,...
                            cell2mat(values(com_map.(fs).(pref).com,curr_key.')));
                    end
                    for ff=["s1d3","s1d6","s2d3","s2d6"]
                        imdata.(pref).(ff)=cat(1,imdata.(pref).(ff),cell2mat(com_map.(fs).(pref).(ff).values(curr_key.')));
                    end
                end

                ffs=["s1d3","s1d6","s2d3","s2d6"];
                for ffidx=1:4
                    if startsWith(pref,"olf") && endsWith(ffs(ffidx),'d6')
%                         keyboard();
                        [~,sortidx]=sort(imdata.(pref).com6);
                    else
                        [~,sortidx]=sort(imdata.(pref).com);
                    end
                    nexttile((prefidx-1)*4+ffidx);
                    if opt.xlim==3
                        plotOne(imdata.(pref).(ffs(ffidx))(sortidx,1:12),'xlim',3,'scale',opt.scale,'gauss2d',opt.gauss2d);
                    else
                        plotOne(imdata.(pref).(ffs(ffidx))(sortidx,:),'scale',opt.scale,'gauss2d',opt.gauss2d);
                    end
                end
            end
        case 3
            fh.("wave"+string(plot_id))=figure('Color','w','Position',[32,32,800,400]);
            tiledlayout(2,2)
            tags=["olf_s1","olf_s2"];
            for prefidx=1:2
                pref=tags(prefidx);
                imdata.(pref)=struct();
                imdata.(pref).com=[];
                for ff=["s1","s2" ]
                    imdata.(pref).(ff)=[];
                end
                for fs1=fss
                    fs=char(fs1);
                    if ~isfield(com_map.(fs),pref)
                        continue
                    end
                    curr_key=com_map.(fs).(pref).com.keys();
                    imdata.(pref).com=cat(1,imdata.(pref).com,...
                        cell2mat(values(com_map.(fs).(pref).com,curr_key.')));
                    for ff=["s1","s2"]
                        imdata.(pref).(ff)=cat(1,imdata.(pref).(ff),cell2mat(com_map.(fs).(pref).(ff).values(curr_key.')));
                    end
                end

                ffs=["s1","s2"];
                for ffidx=1:2
                    [~,sortidx]=sort(imdata.(pref).com);
                    nexttile((prefidx-1)*2+ffidx);
                    plotOne(imdata.(pref).(ffs(ffidx))(sortidx,:),'scale',opt.scale,'gauss2d',opt.gauss2d);
                end
            end
            
    end
end
end

function plotOne(imdata,opt)
arguments
    imdata (:,:) double
    opt.scale (1,2) double = [-1,1]
    opt.title (1,:) char = []
    opt.cmap (1,:) char = 'turbo'
    opt.xlim (1,1) double {mustBeMember(opt.xlim,[3 6])} = 6
    opt.gauss2d (1,1) logical = false
end

colormap(opt.cmap);
if opt.gauss2d
    gk = fspecial('gaussian', [9 3], 1);
    imagesc(conv2(imdata,gk,'same'),opt.scale)
else
    imagesc(imdata,opt.scale)
end
if size(imdata,2)>20
    set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
    xline(12.5,'--w','LineWidth',1);
    xlim([0.5,size(imdata,2)+0.5])
else
    if opt.xlim==6
        xlim([0.5,size(imdata,2)*2+0.5])
    else
        xlim([0.5,size(imdata,2)+0.5])
    end
    set(gca(),'XTick',[0.5,12.5],'XTickLabel',[0,3]);
end
% colorbar();
ylim([0.5,size(imdata,1)+0.5])
if numel(opt.title)>0
    title(opt.title)
end
set(gca(),'YDir','normal')
end