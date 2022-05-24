function fh=sens_dur_wave_bar(sens_map_cells,dur_map_cells,fcom)
% keyboard()
%%

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

diffs=cell(0);
for src=reshape(src_ccfid,1,[])
    t_map=containers.Map( ...
        cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
        sink_src_mat(:,src_ccfid==src));
    sr=map_match(sens_density_map,t_map);
    dr=map_match(dur_density_map,t_map);
    diffs=[diffs;src,idmap.ccfid2reg(src),num2cell(sr),num2cell(dr),num2cell(sr-dr)];
end
[~,idces]=sort(abs(cell2mat(diffs(:,5))),'descend');
diffs=diffs(idces,:);

fh=figure('Color','w','Position',[2,180,900,300]);
tiledlayout(1,3)

[rd,pd]=map_match(dur_density_map,sens_density_map);

for ii=1:5
    r(ii*2-1)=diffs{ii,3};
    r(ii*2)=diffs{ii,4};
end

% for ii=1:5
%     r(ii*2+9)=diffs{end-ii+1,3};
%     r(ii*2+10)=diffs{end-ii+1,4};
% end

nexttile();
hold on
bhd=bar(1,rd,'FaceColor',[0.5,0.5,0.5]);
rs=reshape(r,2,[]).';
bh=bar(2:6,rs,'grouped');
bh(1).FaceColor='w';
bh(2).FaceColor='k';
legend(bh,{'Sensory','Duration'},'Location','northoutside','Orientation','horizontal');
% set(gca(),'XTick',1:10,'XTickLabel',[diffs(1:5,2).',diffs(end-4:end,2).'],'XTickLabelRotation',90)
set(gca(),'XTick',1:6,'XTickLabel',['sens_dur',diffs(1:5,2).'],'XTickLabelRotation',90)
ylabel('Correlation (Spearman''s r)');
ylim([-1,1]);
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
        sr=map_match(sens_TCOM_map,t_map);
        dr=map_match(dur_TCOM_map,t_map);
        diffs=[diffs;src,idmap.ccfid2reg(src),num2cell(sr),num2cell(dr),num2cell(sr-dr)];
    end
    [~,idces]=sort(abs(cell2mat(diffs(:,5))),'descend');
    diffs=diffs(idces,:);

    % sensTCOM, durTCOM vs:  % data from fcom
    % sensDensity, durDensity, connectivity
    % 6s

    [r(1),p(1)]=map_match(sens_TCOM_map,sens_density_map);
    [r(2),p(2)]=map_match(dur_TCOM_map,sens_density_map);
    [r(3),p(3)]=map_match(sens_TCOM_map,dur_density_map);
    [r(4),p(4)]=map_match(dur_TCOM_map,dur_density_map);

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
    bh=bar(rs,'grouped');
    bh(1).FaceColor='w';
    bh(2).FaceColor='k';
    legend(bh,{'Sensory','Duration'},'Location','northoutside','Orientation','horizontal');
%     set(gca(),'XTick',1:12,'XTickLabel',[{'Sensory density','Duration density'},diffs(1:5,2).',diffs(end-4:end,2).'],'XTickLabelRotation',90)
    set(gca(),'XTick',1:7,'XTickLabel',[{'Sensory density','Duration density'},diffs(1:5,2).'],'XTickLabelRotation',90)
    ylabel('Correlation (Spearman''s r)');
    title(['TCOM ',char(tcom_delay)])
    ylim([-1,1]);
    % exportgraphics(gcf(),'collections.pdf','ContentType','vector','Append',true)
end
end
function [r,p]=map_match(map1,map2)
keys=intersect(map1.keys,map2.keys);
v1=cell2mat(values(map1,keys));
v2=cell2mat(values(map2,keys));
[r,p]=corr(v1(:),v2(:),'type','Spearman');
end