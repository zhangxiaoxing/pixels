% Borrowed from wave.sens_dur_TCOM_corr, consider merge
function [fh,out]=reg_wave_movie(fcom,opt)
arguments
    fcom
    opt.spatial_coord (1,1) logical = true
    opt.write_video (1,1) logical = false
    opt.export_data (1,1) logical = false
    opt.color_type (1,:) char {mustBeMember(opt.color_type,{'proportion','wavefront'})} = 'proportion'
    opt.annotation (1,:) char {mustBeMember(opt.annotation,{'500ms','std'})} = 'std'
    opt.slow_rate (1,1) double {mustBeMember(opt.slow_rate,[4,8,10])} = 8
    opt.waves (1,:) char {mustBeMember(opt.waves,{'sens','dur','both'})} = both
end
if opt.spatial_coord
    coordmap=bz.gephi.reg2coord();
end
if opt.spatial_coord
    fnsuffix="_spatial_"+string(opt.waves);
else
    fnsuffix="_"+string(opt.waves);
end
if opt.write_video
    v=VideoWriter("sens_dur_wave"+fnsuffix+".mp4",'MPEG-4');
    open(v);
end
fh=figure('Color','w','Position',[16,180,740,480]);

% dxmap=containers.Map({'AON'},[0.02])



% [rs,rp]=deal([]);
for dd=["d3"] %["d3","d6"]
    inter_reg=intersect(...
        intersect(ephys.getGreyRegs('range','grey'),...
        fcom.(dd).collection(:,2)),...
        fcom.dur.collection(:,2));

    [~,sidx]=ismember(inter_reg,fcom.(dd).collection(:,2));
    [~,didx]=ismember(inter_reg,fcom.dur.collection(:,2));

    sv=cell2mat(fcom.(dd).collection(sidx,[1,5]));
    dv=cell2mat(fcom.dur.collection(didx,[1,5]));
    
    [~,yidx]=sort(dv(:,1));
    [~,xidx]=sort(sv(:,1));

%     coord=nan(numel(inter_reg),3);
    regs=cell(numel(inter_reg),1);
    out=cell(numel(inter_reg)+1,332);
    out(2:end,1)=inter_reg;
    out{1,1}='time';

    tidx=2;
    switch opt.slow_rate
        case 4
            dT=1/30;
        case 8
            dT=1/60;
        case 10
            dT=1/75;
    end

    for timebin=dT:(dT):12
        disp(timebin)
        clf
        hold on;
        if ismember(opt.waves,{'both','sens'})
            xline(timebin/4,'r-')
        end
        if ismember(opt.waves,{'both','dur'})
            yline(timebin/4,'b-')
        end
        %TODO std envelope
        switch opt.annotation
            case '500ms'
                rfh=fill((timebin+[-1,1,1,-1])./4,[0,0,3,3],'r','EdgeColor','none','FaceAlpha',0.1);
                bfh=fill([0,0,3,3],(timebin+[-1,1,1,-1])./4,'b','EdgeColor','none','FaceAlpha',0.1);
            case 'std'
% 
%                 int_dv=interp1([0;dv(yidx,1);12],[sv(yidx(1),2);sv(yidx,2);sv(yidx(end),2)],0:0.1:12);
%                 int_sv=interp1([0;sv(xidx,1);12],[dv(xidx(1),2);dv(xidx,2);dv(xidx(end),2)],0:0.1:12);
%                 fill((timebin+[-0.5*int_sv,0.5*fliplr(int_dv)])./4,[0:0.1:12,12:-0.1:0]./4,'r','EdgeColor','none','FaceAlpha',0.1);
%                 fill([0:0.1:12,12:-0.1:0]./4,(timebin+[-0.5*int_dv,0.5*fliplr(int_dv)])./4,'b','EdgeColor','none','FaceAlpha',0.1);

                svs=[sv(yidx(1),2);sv(yidx,2);sv(yidx(end),2)];
                sxs=[0;dv(yidx,1);12];
                dvs=[dv(xidx(1),2);dv(xidx,2);dv(xidx(end),2)];
                dxs=[0;sv(xidx,1);12];
                if ismember(opt.waves,{'both','sens'})
                    rfh=fill((timebin+[-0.5*svs;0.5*flip(svs)])./4,[sxs;flip(sxs)]./4,'r','EdgeColor','none','FaceAlpha',0.1);
                end
                if ismember(opt.waves,{'both','dur'})    
                    bfh=fill([dxs;flip(dxs)]./4,(timebin+[-0.5*dvs;0.5*flip(dvs)])./4,'b','EdgeColor','none','FaceAlpha',0.1);
                end

        end

        for ri=1:numel(inter_reg)
            sens_idx=(strcmp(fcom.(dd).collection(:,2),inter_reg(ri)));
            dur_idx=(strcmp(fcom.dur.collection(:,2),inter_reg(ri)));
            if opt.spatial_coord
                coord=coordmap(inter_reg{ri});
                xx=coord(1);
                yy=coord(2);
                txx=xx;
            else
                yy=fcom.dur.collection{dur_idx,1}./4;
                xx=fcom.(dd).collection{sens_idx,1}./4;
                txx=xx;
                if ismember(inter_reg(ri),{'DP','RSP','PL'})
                    txx=txx-0.015;
                elseif ismember(inter_reg(ri),{'SS','ILA'})
                    txx=txx-0.025;
                elseif ismember(inter_reg(ri),{'ORB','GPe','VIS','LAT'})
                    txx=txx+0.015;
                end
                
            end
%             coord(ri,:)=[xx,yy,1];
            regs(ri)=inter_reg(ri);
            switch opt.color_type
                case 'proportion'
                    sens_sel=strcmp(fcom.(dd).com_meta(:,8),inter_reg(ri));
                    sens_TCOM=cell2mat(fcom.(dd).com_meta(sens_sel,3));
                    sens_prop=nnz(sens_TCOM<=timebin & sens_TCOM>(timebin-1))./numel(sens_TCOM);

                    dur_sel=strcmp(fcom.dur.com_meta(:,8),inter_reg(ri));
                    dur_TCOM=cell2mat(fcom.dur.com_meta(dur_sel,3));
                    dur_prop=nnz(dur_TCOM<=timebin & dur_TCOM>(timebin-1))./numel(dur_TCOM);

                    sens_c=min([(sens_prop./0.2).*0.8,0.8]);
                    dur_c=min([(dur_prop./0.2).*0.8,0.8]);
                case 'wavefront'
                    sens_sel=strcmp(fcom.(dd).collection(:,2),inter_reg(ri));
                    sens_TCOM=cell2mat(fcom.(dd).collection(sens_sel,[1,5]));
                    if strcmp(opt.annotation,'500ms')
                        if (sens_TCOM(1)-1)<=timebin && (sens_TCOM(1)+1)>timebin
                            sens_c=0.8;
                        else
                            sens_c=0;
                        end
    
                        dur_sel=strcmp(fcom.dur.collection(:,2),inter_reg(ri));
                        dur_TCOM=cell2mat(fcom.dur.collection(dur_sel,[1,5]));
                        if (dur_TCOM(1)-1)<=timebin && (dur_TCOM(1)+1)>timebin
                            dur_c=0.8;
                        else
                            dur_c=0;
                        end
                    elseif strcmp(opt.annotation,'std')
                        if (sens_TCOM(1)-0.5*sens_TCOM(2))<=timebin && (sens_TCOM(1)+0.5*sens_TCOM(2))>timebin
                            sens_c=0.8;
                        else
                            sens_c=0;
                        end

                        dur_sel=strcmp(fcom.dur.collection(:,2),inter_reg(ri));
                        dur_TCOM=cell2mat(fcom.dur.collection(dur_sel,[1,5]));
                        if (dur_TCOM(1)-0.5*dur_TCOM(2))<=timebin && (dur_TCOM(1)+0.5*dur_TCOM(2))>timebin
                            dur_c=0.8;
                        else
                            dur_c=0;
                        end

                    end
                    
            end
            switch opt.waves
                case 'both'
                    wave_color=[0.2+sens_c,0.2,0.2+dur_c];
                case 'sens'
                    wave_color=[0.2+sens_c,0.2,0.2];
                case 'dur'
                    wave_color=[0.2,0.2,0.2+dur_c];
            end

            scatter(xx,yy,49,'o','MarkerFaceColor',wave_color,'MarkerEdgeColor','none');
            text(txx,yy-0.002,inter_reg{ri},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',12,'Color',wave_color);
            out{1,tidx}=timebin;
            out{ri+1,tidx}=[sens_c>0,dur_c>0];
        end
        tidx=tidx+1;
        ah0=gca();
        title(sprintf('%0.2fs',timebin/4),'FontSize',12);
        set(ah0,'TitleHorizontalAlignment','right')

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
        switch opt.waves
            case 'both'
                legend([rfh,bfh],{'Sensory wave','Duration wave'},'Location','northeastoutside')
            case 'sens'
                legend([rfh],{'Sensory wave'},'Location','northeastoutside')
            case 'dur'
                legend([bfh],{'Duration wave'},'Location','northeastoutside')
        end
        set(ah0,'Position',[0.1,0.15,0.6,0.7])
        if opt.write_video
            writeVideo(v,getframe(fh))
        end
    end
    if opt.write_video
        close(v)
    end
end

if opt.export_data
    if exist('wave_movie.hdf5','file')
        delete('wave_movie.hdf5')
    end
    wave_v=cell2mat(out(2:end,2:end));
    sens_v=wave_v(:,1:2:end);
    dur_v=wave_v(:,2:2:end);
    reg=out(2:end,1);
    h5create('wave_movie.hdf5','/sens_v',size(sens_v),'Datatype','double')
    h5write('wave_movie.hdf5','/sens_v',double(sens_v))
    h5create('wave_movie.hdf5','/dur_v',size(dur_v),'Datatype','double')
    h5write('wave_movie.hdf5','/dur_v',double(dur_v))
    h5create('wave_movie.hdf5','/reg',size(reg),'Datatype','string')
    h5write('wave_movie.hdf5','/reg',reg)
end
end