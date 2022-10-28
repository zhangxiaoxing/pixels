
function fh=plot_pct_wave(pct_meta,com_map,opt)
arguments
%     eff_meta
    pct_meta
    com_map
    opt.comb_set (1,:) double {mustBeInteger,mustBePositive} = 1
    opt.sort_by (1,1) double {mustBeMember(opt.sort_by,[3 6])} = 6
    opt.xlim (1,1) double {mustBeMember(opt.xlim,[3 6])} = 6
    opt.lesser_grade (1,1) logical = false
    opt.scale (1,2) double = [-1,1]

end
%% global

%     %TODO do not repeat>>>>>>>>>>>>>>>>
%     sens_efsz=max(abs(eff_meta.cohen_d_olf),[],2);
%     sens_win=[min(sens_efsz)./2,prctile(sens_efsz,[20:20:100])];
% 
%     dur_efsz=max(abs(eff_meta.cohen_d_dur),[],2);
%     dur_win=[min(dur_efsz)./2,prctile(dur_efsz,[20:20:100])];
%     % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% 
%     pct_meta=pct.get_pct_meta(eff_meta,sens_efsz,sens_win,dur_efsz,dur_win);
% 
%     com_map=wave.get_pct_com_map(pct_meta,'curve',true);


for plot_id=opt.comb_set
    fss=reshape(fieldnames(com_map),1,[]);
    imdata=struct();
    switch plot_id
        case 1 % multiplexed
            fh.("wave"+string(plot_id))=figure('Color','w','Position',[32,32,1400,800]);
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
                    if opt.sort_by==6
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
                        plotOne(imdata.(pref).(tags(ffidx))(sortidx,1:12),'xlim',3,'scale',opt.scale);
                    else
                        plotOne(imdata.(pref).(tags(ffidx))(sortidx,:),'scale',opt.scale);
                    end
                end
            end
        case 2 % olf
            fh.("wave"+string(plot_id))=figure('Color','w','Position',[32,32,1400,800]);
            tiledlayout(4,4)
            tags=["olf_s1","olf_s2","dur_d3","dur_d6"];
            for prefidx=1:4
                pref=tags(prefidx);
                imdata.(pref)=struct();
                imdata.(pref).com=[];
                for ff=["s1d3","s1d6","s2d3","s2d6"]
                    imdata.(pref).(ff)=[];
                end
                for fs1=fss
                    fs=char(fs1);
                    if opt.sort_by==6
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
                ffs=["s1d3","s1d6","s2d3","s2d6"];
                for ffidx=1:4
                    nexttile((prefidx-1)*4+ffidx);
                    if opt.xlim==3
                        plotOne(imdata.(pref).(ffs(ffidx))(sortidx,1:12),'xlim',3,'scale',opt.scale);
                    else
                        plotOne(imdata.(pref).(ffs(ffidx))(sortidx,:),'scale',opt.scale);
                    end
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
    opt.cmap (1,:) char = 'bone'
    opt.xlim (1,1) double {mustBeMember(opt.xlim,[3 6])} = 6
end

colormap(opt.cmap);
if false
    gk = fspecial('gaussian', [3 3], 1);
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
colorbar();
ylim([0.5,size(imdata,1)+0.5])
if numel(opt.title)>0
    title(opt.title)
end
set(gca(),'YDir','normal')
end