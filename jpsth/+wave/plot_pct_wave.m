
function fh=plot_pct_wave(eff_meta,opt)
arguments
    eff_meta
    opt.comb_set (1,:) double {mustBeInteger,mustBePositive} = 1
end
%% global
for plot_id=opt.comb_set
    %TODO do not repeat>>>>>>>>>>>>>>>>
    sens_efsz=max(abs(eff_meta.cohen_d_olf),[],2);
    sens_win=[min(sens_efsz)./2,prctile(sens_efsz,[20:20:100])];

    dur_efsz=max(abs(eff_meta.cohen_d_dur),[],2);
    dur_win=[min(dur_efsz)./2,prctile(dur_efsz,[20:20:100])];
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    pct_meta4=pct.get_pct_meta(eff_meta,sens_efsz,sens_win,dur_efsz,dur_win,'single_mod_thresh',4);
    com_map=wave.get_pct_com_map(pct_meta4,'curve',true);
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
                    curr_key=com_map.(fs).(pref).com.keys();
                    imdata.(pref).com=cat(1,imdata.(pref).com,cell2mat(values(com_map.(fs).(pref).com,curr_key.')));
                    for ff=["s1d3","s1d6","s2d3","s2d6"]
                        imdata.(pref).(ff)=cat(1,imdata.(pref).(ff),cell2mat(com_map.(fs).(pref).(ff).values(curr_key.')));
                    end
                end

                [~,sortidx]=sort(imdata.(pref).com);
                for ffidx=1:4
                    nexttile((prefidx-1)*4+ffidx);
                    plotOne(imdata.(pref).(tags(ffidx))(sortidx,:));
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
                    curr_key=com_map.(fs).(pref).com.keys();
                    imdata.(pref).com=cat(1,imdata.(pref).com,cell2mat(values(com_map.(fs).(pref).com,curr_key.')));
                    for ff=["s1d3","s1d6","s2d3","s2d6"]
                        imdata.(pref).(ff)=cat(1,imdata.(pref).(ff),cell2mat(com_map.(fs).(pref).(ff).values(curr_key.')));
                    end
                end

                [~,sortidx]=sort(imdata.(pref).com);
                ffs=["s1d3","s1d6","s2d3","s2d6"];
                for ffidx=1:4
                    nexttile((prefidx-1)*4+ffidx);
                    plotOne(imdata.(pref).(ffs(ffidx))(sortidx,:));
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
end

colormap(opt.cmap);
gk = fspecial('gaussian', [3 3], 1);
imagesc(conv2(imdata,gk,'same'),opt.scale)

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
set(gca(),'YDir','normal')
end