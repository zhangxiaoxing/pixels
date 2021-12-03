function rtn=plot_global_wave_sorted(opt)
arguments
    opt.plot_global_wave (1,1) logical = false
    opt.plot_sel_idx (1,1) logical = true
end

if opt.plot_global_wave
    immata=[];
    immatb=[];
    antiTCOM=[];
    TCOM=[];
    for sutype=["any_s1", "any_s2","any_nonmem"]
        com_map_6=wave.get_com_map('curve',true,'rnd_half',false,'delay',6,'cell_type',sutype);
        for fsc=reshape(fieldnames(com_map_6),1,[])
            fs=fsc{1};
            samp_key_S1=cell2mat(com_map_6.(fs).s1.keys);
            samp_key_S2=cell2mat(com_map_6.(fs).s2.keys);
            COMS1=cell2mat(values(com_map_6.(fs).s1,num2cell(samp_key_S1)));
            COMS2=cell2mat(values(com_map_6.(fs).s2,num2cell(samp_key_S2)));

            sortmat=[ones(size(COMS1)),2*ones(size(COMS2));...
                double(samp_key_S1),double(samp_key_S2);...
                COMS1,COMS2].';

            for ri=1:size(sortmat,1)
                one_heat=com_map_6.(fs).(sprintf('s%dcurve',sortmat(ri,1)))(sortmat(ri,2));
                immata=[immata;one_heat];
                TCOM=[TCOM;com_map_6.(fs).(sprintf('s%d',sortmat(ri,1)))(sortmat(ri,2))];
            end

            for ri=1:size(sortmat,1)
                one_heat=com_map_6.(fs).(sprintf('s%danticurve',sortmat(ri,1)))(sortmat(ri,2));
                immatb=[immatb;one_heat];
            end
        end


        normmax=max(abs([immata,immatb]),[],2);
        immatan=immata./normmax;
        immatbn=immatb./normmax;
        [TCOM,tcomidx]=sort(TCOM);
        %% actual plot
        fh=figure('Color','w','Position',[32,32,750,350]);
        plotOne(1,immatan(tcomidx,:));
        plotOne(2,immatbn(tcomidx,:));
        exportgraphics(fh,sprintf('global_wave_%s.pdf',sutype),'ContentType','vector')
    end
    rtn=TCOM
end

%% global selectivity
if opt.plot_sel_idx
    com_map=wave.get_com_map('curve',true,'rnd_half',false,'delay',6,'cell_type','ctx_trans','selidx',true);

    sortmat_all=[];
    for fsc=reshape(fieldnames(com_map),1,[])
        fs=fsc{1};
        samp_key_S1=cell2mat(com_map.(fs).s1.keys);
        samp_key_S2=cell2mat(com_map.(fs).s2.keys);
        COMS1=cell2mat(values(com_map.(fs).s1,num2cell(samp_key_S1)));
        COMS2=cell2mat(values(com_map.(fs).s2,num2cell(samp_key_S2)));

        sortmat=[ones(size(COMS1)),2*ones(size(COMS2));...
            double(samp_key_S1),double(samp_key_S2);...
            COMS1,COMS2].';

        immat=[];
        for ri=1:size(sortmat,1)
            one_heat=com_map.(fs).(sprintf('s%dcurve',sortmat(ri,1)))(sortmat(ri,2));
            if ~all(isfinite(one_heat),'all')
                keyboard()
            end
            immat=[immat;one_heat];
        end
        sortmat_all=[sortmat_all;sortmat,immat];
    end

    sortmat_all=sortrows(sortmat_all,3);
    rtn=sortmat_all;
    %% actual plot
    fh=figure('Color','w','Position',[32,32,750,350]);
    plotOne(1,sortmat_all(:,4:end),'cmap','copper','scale',[-0.1,0.4],'pnl_num',1);
    xlim([8.5,44.5])
    exportgraphics(fh,'ctx_trans_sel_idx_delay.pdf','ContentType','vector')
end

end


function plotOne(subidx,imdata,opt)
arguments
    subidx (1,1) double {mustBeInteger,mustBePositive}
    imdata double
    opt.cmap (1,:) char = 'turbo'
    opt.scale (1,2) double = [-0.6,0.6]
    opt.pnl_num (1,1) double {mustBeInteger,mustBePositive} = 2
end
subplot(1,opt.pnl_num,subidx);
hold on
if strcmp(opt.cmap,'bluewhitered')
    colormap(bluewhitered(256))
    colormap(bluewhitered(256))
else
    colormap(opt.cmap);
end
gk = fspecial('gaussian', [3 3], 1);
imAlpha=ones(size(imdata));
imAlpha(isnan(imdata))=0;
im=imagesc(conv2(imdata,gk,'same'),opt.scale);
im.AlphaData=imAlpha;
set(gca(),'color',[1 1 1],'XTick',[0.5,20.5,40.5]+16,'XTickLabel',[0,5,10]);
arrayfun(@(x) xline(x,'--w'),[12,16,40,44]+0.5)
colorbar();
ylim([0.5,size(imdata,1)+0.5])
xlim([0.5,size(imdata,2)+0.5])
end