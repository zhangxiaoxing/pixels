function fh=sens_dur_wave_bar(sens_map_cells,dur_map_cells,fcom)
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');

rsel=@(x) cell2mat(x(:,3))==5 & ismember(x(:,2),ephys.getGreyRegs);

sens_density_map=containers.Map( ...
    sens_map_cells{3}.keys(), ...
    cellfun(@(x) x(1), sens_map_cells{3}.values));

dur_density_map=containers.Map( ...
    dur_map_cells{3}.keys(), ...
    cellfun(@(x) x(1), dur_map_cells{3}.values));

% sensDensity, durDensity vs connectivity:
fh=figure('Color','w','Position',[2,2,2400,1200]);
th=tiledlayout(2,6);
for cycle=1:2
    if cycle==1
        loglogtype='Spearman';
        linearlogtype='Spearman';
        ylbl='Spearman''s r';
    else
        loglogtype='PearsonLogLog';
        linearlogtype='PearsonLinearLog';
        ylbl='Pearson''s r';
    end
    diffs=cell(0);
    for src=reshape(src_ccfid,1,[])
        t_map=containers.Map( ...
            cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
            sink_src_mat(:,src_ccfid==src));
        sr=map_match(sens_density_map,t_map,loglogtype);%log log
        dr=map_match(dur_density_map,t_map,loglogtype);%log log
        diffs=[diffs;src,idmap.ccfid2reg(src),num2cell(sr),num2cell(dr),num2cell(sr-dr)];
    end
    [~,idces]=sort(abs(cell2mat(diffs(:,5))),'descend');
    diffs=diffs(idces,:);
    nexttile();
    hold on
    [rd,p]=map_match(dur_density_map,sens_density_map,loglogtype);%log log
    bhd=bar(1,rd,'FaceColor',[0.5,0.5,0.5]);
    set(gca(),'XTick',1,'XTickLabel',{'sens_dur'},'XTickLabelRotation',90,...
        'YLim',[-1,1],'YTick',-1:0.5:1)
    ylabel(ylbl);
    nexttile();
    r=nan(5,1);
    for ii=1:5
        r(ii*2-1)=diffs{ii,3};
        r(ii*2)=diffs{ii,4};
    end


    % for ii=1:5
    %     r(ii*2+9)=diffs{end-ii+1,3};
    %     r(ii*2+10)=diffs{end-ii+1,4};
    % end

    rs=reshape(r,2,[]).';
    bh=bar(1:5,rs,'grouped');
    bh(1).FaceColor='w';
    bh(2).FaceColor='k';
    legend(bh,{'Sensory','Duration'},'Location','northoutside','Orientation','horizontal');
    % set(gca(),'XTick',1:10,'XTickLabel',[diffs(1:5,2).',diffs(end-4:end,2).'],'XTickLabelRotation',90)
    set(gca(),'XTick',1:5,'XTickLabel',diffs(1:5,2).','XTickLabelRotation',90, ...
        'YLim',[-1,1],'YTick',-1:0.5:1)
    ylabel(ylbl);
    title("density")





    for tcom_delay=["d3","d6"]
        sens_TCOM_map=containers.Map( ...
            fcom.(tcom_delay).collection(rsel(fcom.(tcom_delay).collection),2), ...
            fcom.(tcom_delay).collection(rsel(fcom.(tcom_delay).collection),1) ...
            );
        dur_TCOM_map=containers.Map( ...
            fcom.dur.collection(rsel(fcom.dur.collection),2), ...
            fcom.dur.collection(rsel(fcom.dur.collection),1) ...
            );

        diffs=cell(0);
        for src=reshape(src_ccfid,1,[])
            t_map=containers.Map( ...
                cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
                sink_src_mat(:,src_ccfid==src));
            sr=map_match(sens_TCOM_map,t_map,linearlogtype);%linear-log
            dr=map_match(dur_TCOM_map,t_map,linearlogtype);%linear-log
            diffs=[diffs;src,idmap.ccfid2reg(src),num2cell(sr),num2cell(dr),num2cell(sr-dr)];
        end
        [~,idces]=sort(abs(cell2mat(diffs(:,5))),'descend');
        diffs=diffs(idces,:);

        % sensTCOM, durTCOM vs:  % data from fcom
        % sensDensity, durDensity, connectivity
        % 6s

        [r(1),p(1)]=map_match(sens_TCOM_map,sens_density_map,linearlogtype);%linear-log
        [r(2),p(2)]=map_match(dur_TCOM_map,sens_density_map,linearlogtype);%linear-log
        [r(3),p(3)]=map_match(sens_TCOM_map,dur_density_map,linearlogtype);%linear-log
        [r(4),p(4)]=map_match(dur_TCOM_map,dur_density_map,linearlogtype);%linear-log

        for ii=1:5
            r(ii*2+3)=diffs{ii,3};
            r(ii*2+4)=diffs{ii,4};
        end

        %     for ii=1:5
        %         r(ii*2+13)=diffs{end-ii+1,3};
        %         r(ii*2+14)=diffs{end-ii+1,4};
        %     end
        nexttile();
        rs=reshape(r,2,[]).';
        bh=bar(rs(1:2,:),'grouped');
        bh(1).FaceColor='w';
        bh(2).FaceColor='k';
        set(gca(),'XTick',1:2,'XTickLabel',{'Sensory density','Duration density'},'XTickLabelRotation',90)
        ylim([-1,1]);
        ylabel(ylbl);
        nexttile()
        bh=bar(rs(3:end,:),'grouped');
        bh(1).FaceColor='w';
        bh(2).FaceColor='k';
        legend(bh,{'Sensory','Duration'},'Location','northoutside','Orientation','horizontal');
        %     set(gca(),'XTick',1:12,'XTickLabel',[{'Sensory density','Duration density'},diffs(1:5,2).',diffs(end-4:end,2).'],'XTickLabelRotation',90)
        set(gca(),'XTick',1:5,'XTickLabel',diffs(1:5,2).','XTickLabelRotation',90,...
            'YLim',[-1,1],'YTick',-1:0.5:1)
        ylabel(ylbl);
        title(['TCOM ',char(tcom_delay)])
    end
    % exportgraphics(gcf(),'collections.pdf','ContentType','vector','Append',true)
end
end
function [r,p]=map_match(map1,map2,corr_type)
arguments
    map1
    map2
    corr_type (1,:) char {mustBeMember(corr_type,{'Spearman','PearsonLogLog','PearsonLinearLog'})}
end
keys=intersect(map1.keys,map2.keys);
v1=cell2mat(values(map1,keys));
v2=cell2mat(values(map2,keys));
switch corr_type
    case 'Spearman'
        [r,p]=corr(v1(:),v2(:),'type','Spearman');
    case 'PearsonLogLog'
        vsel=v1>0 & v2>0;
        [r,p]=corr(log10(v1(vsel)).',log10(v2(vsel)).','type','Pearson');
    case 'PearsonLinearLog'
        vsel=v2>0;
        [r,p]=corr(v1(vsel).',log10(v2(vsel)).','type','Pearson');
end

end