% Borrowed from wave.sens_dur_TCOM_corr, consider merge
function [fh,out]=reg_wave_movie(fcom,opt)
arguments
    fcom
    opt.spatial_coord (1,1) logical = true
end
if opt.spatial_coord
    coordmap=bz.gephi.reg2coord();
end
if opt.spatial_coord
    fnsuffix="_spatial";
else
    fnsuffix="";
end

v=VideoWriter("sens_dur_wave"+fnsuffix+".mp4",'MPEG-4');
open(v);
fh=figure('Color','w','Position',[16,180,1280,720]);

% dxmap=containers.Map({'AON'},[0.02])




% [rs,rp]=deal([]);
for dd=["d3"] %["d3","d6"]
    inter_reg=intersect(...
        intersect(ephys.getGreyRegs('range','grey'),...
        fcom.(dd).collection(:,2)),...
        fcom.dur.collection(:,2));

%     coord=nan(numel(inter_reg),3);
    regs=cell(numel(inter_reg),1);
    out=cell(numel(inter_reg)+1,332);
    out(2:end,1)=inter_reg;
    out{1,1}='time';

    tidx=2;
    for timebin=1:(1/30):12
        disp(timebin)
%         keyboard()
%         pause(0.1)
        clf
        %TODO insert video frame handling
        hold on;

        for ri=1:numel(inter_reg)
            sens_idx=(strcmp(fcom.(dd).collection(:,2),inter_reg(ri)));
            dur_idx=(strcmp(fcom.dur.collection(:,2),inter_reg(ri)));
            if opt.spatial_coord
                coord=coordmap(inter_reg{ri});
                xx=coord(1);
                yy=coord(2);
            else
                yy=fcom.dur.collection{dur_idx,1}./4;
                xx=fcom.(dd).collection{sens_idx,1}./4;
            end
%             coord(ri,:)=[xx,yy,1];
            regs(ri)=inter_reg(ri);
            %TODO dynamically adjust color
            sens_sel=strcmp(fcom.(dd).com_meta(:,8),inter_reg(ri));
            sens_TCOM=cell2mat(fcom.(dd).com_meta(sens_sel,3));
            sens_prop=nnz(sens_TCOM<=timebin & sens_TCOM>(timebin-1))./numel(sens_TCOM);

            dur_sel=strcmp(fcom.dur.com_meta(:,8),inter_reg(ri));
            dur_TCOM=cell2mat(fcom.dur.com_meta(dur_sel,3));
            dur_prop=nnz(dur_TCOM<=timebin & dur_TCOM>(timebin-1))./numel(dur_TCOM);

            sens_c=min([(sens_prop./0.2).*0.8,0.8]);
            dur_c=min([(dur_prop./0.2).*0.8,0.8]);

            wave_color=[0.2+sens_c,0.2,0.2+dur_c];

            scatter(xx,yy,64,'o','MarkerFaceColor',wave_color,'MarkerEdgeColor','none');
            text(xx,yy-0.005,inter_reg{ri},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',16,'Color',wave_color);
            out{1,tidx}=timebin;
            out{ri+1,tidx}=[sens_prop,dur_prop];
            tidx=tidx+1;

        end
        ah0=gca();
        title(sprintf('%0.2fs',timebin/4),'FontSize',12);
        set(ah0,'TitleHorizontalAlignment','right')
%         colorbar(gca(),"eastoutside",)

        ax1=gca();
        ax2=axes();
        linkaxes([ax1,ax2]);
        ax2.Visible='off';
        [ax2.XTick,ax2.YTick]=deal([]);
        redmap=zeros(64,3)+0.2+(0:(0.8/63):0.8).'*[1,0,0];
        bluemap=zeros(64,3)+0.2+(0:(0.8/63):0.8).'*[0,0,1];
        colormap(ax1,redmap);
        colormap(ax2,bluemap);
        set([ax1,ax2],'Position',[0.1,0.15,0.6,0.7])
        ah1=colorbar(ax1,'Position',[0.75,0.15,0.03,0.7],'Limits',[0,1]);
        ah2=colorbar(ax2,'Position',[0.85,0.15,0.03,0.7],'Limits',[0,1]);
        set([ah1,ah2],'Ticks',[0,0.5,1]);
        set([ah1,ah2],'TickLabels',[0,10,20]);
        ah1.Label.FontSize=12;
        ah2.Label.FontSize=12;
        ah1.Label.String='Proportion of neuron at odor TCOM (%)';
        ah2.Label.String='Proportion of neuron at duration TCOM (%)';
        if opt.spatial_coord
            xlabel(ah0,'AP (\mum)')
            ylabel(ah0,'DV (\mum)')
            ylim(ah0,[75,575])
            xlim(ah0,[250,900])
            set(ah0,'YDir','reverse')
        else
            xlabel(ah0,'Sensory wave timing (s)')
            ylabel(ah0,'Duration wave timing (s)')
            xlim(ah0,[1.25,2.06])
            ylim(ah0,[1.25,2.06])
        end
        writeVideo(v,getframe(fh))
    end
    close(v)
end
end