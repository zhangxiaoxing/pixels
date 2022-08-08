
function fh=plot_pct_wave(eff_meta,opt)
arguments
    eff_meta
    opt.comb_set (1,:) double {mustBeInteger,mustBePositive} = 1
end
%% global
for plot_id=opt.comb_set
    switch plot_id
        case 1 % both
            com_map=wave.get_pct_com_map(eff_meta,'curve',true);
            tag='Mixed selective';
            fn='Mixed_selective_wave_heat.pdf';
%         case 2 % sens_only
%             com_map_3=wave.get_com_map(sel_meta,'curve',true,'rnd_half',false,'delay',3,'wave','onlyContext2');
%             com_map_6=wave.get_com_map(sel_meta,'curve',true,'rnd_half',false,'delay',6,'wave','onlyContext2');%             com_map_6_si=com_map_6_sel_si;
%             tag='Selective only in 6s trials';
%             fn='corr_3_6_only6.pdf';
%         case 3 % dur_only
%             com_map_3=wave.get_com_map(sel_meta,'curve',true,'rnd_half',false,'delay',3,'wave','onlyContext1');
%             com_map_6=wave.get_com_map(sel_meta,'curve',true,'rnd_half',false,'delay',6,'wave','onlyContext1');
%             tag='Selective only in 3s trials';
%             fn='corr_3_6_only3.pdf';
    end
    fss=reshape(fieldnames(com_map),1,[]);
    imdata=struct();
    for pref=["s1d3","s1d6","s2d3","s2d6"]
        imdata.(pref)=struct();
        imdata.(pref).com=[];
        for ff=["s1d3","s1d6","s2d3","s2d6"]
            imdata.(pref).(ff)=[];
        end
    end
%     [im_s1d3,im_s2d3,im_s1d6,im_s2d6]=deal([]);
    
    for fs1=fss
        fs=char(fs1);
        for pref=["s1d3","s1d6","s2d3","s2d6"]
            curr_key=com_map.(fs).(pref).com.keys();
            imdata.(pref).com=cat(1,imdata.(pref).com,cell2mat(values(com_map.(fs).(pref).com,curr_key.')));
            for ff=["s1d3","s1d6","s2d3","s2d6"]
                imdata.(pref).(ff)=cat(1,imdata.(pref).(ff),cell2mat(com_map.(fs).(pref).(ff).values(curr_key.')));
            end
        end
    end

    %% wave
    fh.("wave"+string(plot_id))=figure('Color','w','Position',[32,32,1400,800]);
    tiledlayout(4,4)
    tags=["s1d3","s1d6","s2d3","s2d6"];
    for prefidx=1:4
        [~,sortidx]=sort(imdata.(tags(prefidx)).com);
        for ffidx=1:4
            nexttile((prefidx-1)*4+ffidx);
            plotOne(imdata.(tags(prefidx)).(tags(ffidx))(sortidx,:));
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