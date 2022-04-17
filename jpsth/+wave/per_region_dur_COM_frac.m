function [fcom,ffrac]=per_region_dur_COM_frac(opt)
arguments
    opt.frac_COM (1,1) logical = true
    opt.frac_PVSST (1,1) logical = false
    opt.COM_PVSST (1,1) logical = false
    opt.frac_sensemotor (1,1) logical = false
    opt.COM_sensemotor (1,1) logical = true
    opt.sust_type (1,:) char {mustBeMember(opt.sust_type,{'any','sust','trans'})} = 'any'
    opt.range (1,:) char {mustBeMember(opt.range,{'CH','CTX','grey'})} = 'grey'
    opt.corr (1,:) char {mustBeMember(opt.corr,{'Pearson','Spearman'})} = 'Spearman'
    opt.export (1,1) logical = true
    opt.selidx (1,1) logical = false % calculate COM of selectivity index
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6 % DPA delay duration
end

fn=sprintf('d%d',opt.delay);
ureg=ephys.getGreyRegs('range',opt.range).';
%data of interest, region,branch level, count
com_str=wave.get_dur_com_map();% ffrac.collection=ephys.per_region_fraction('memtype',opt.sust_type,'delay',opt.delay);
meta=ephys.util.load_meta();
tempstats=cell(0,0);
for sess=reshape(fieldnames(com_str),1,[])
    sessid=str2double(replace(sess,'s',''));
    idkey=uint16(cell2mat(com_str.(sess{1}).(fn).keys()));
    metaid=meta.allcid(meta.sess==sessid);
    [~,meta_loc]=ismember(idkey,metaid);
    comreg=subsref(meta.reg_tree(5,meta.sess==sessid),struct(type='()',subs={{meta_loc}})).';
    tempstats=[tempstats;com_str.(sess{1}).(fn).values().',comreg,num2cell(repmat(sessid,numel(idkey),1)),num2cell(idkey).'];
end
fcom.collection=cell(0);
for reg=reshape(ureg,1,[])
    regsel=strcmp(tempstats(:,2),reg);
    if any(regsel)
        fcom.collection=[fcom.collection;{mean([tempstats{regsel,1}]),reg{1},5,nnz(regsel)}];
    end
end

%% fracmap
waveid=ephys.get_wave_id(meta.sess,meta.allcid);
anovameta=wave.get_dur_waveid();
%sel for wave id 7&8
waveid(waveid==0 & anovameta.dur_waveid==3)=7;
waveid(waveid==0 & anovameta.dur_waveid==6)=8;
ffrac.collection=cell(0,0);%cell(numel(regs),4)
dursel=waveid>6;
for reg=reshape(ureg,1,[])
    senssel=nnz(dursel&strcmp(meta.reg_tree(5,:),reg).');
    regsel=nnz(strcmp(meta.reg_tree(5,:),reg).');
    ffrac.collection=[ffrac.collection;{senssel/regsel,reg{1},5,regsel}];
end
%% hiermap
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
    src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
    sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');
    anov=sink_src_mat(:,src_ccfid==idmap.reg2ccfid('LAT'));
    hiermap=containers.Map(cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
        anov);

%% com vs frac
if opt.frac_COM
    fh=figure('Color','w','Position',[100,100,245,235]);
    hold on;
    coord=[];
    regs=[];
    inter_reg=intersect(fcom.collection(:,2),ureg);
    for ri=1:numel(inter_reg)
        fridx=find(strcmp(ffrac.collection(:,2),inter_reg(ri)));
        if ~isempty(fridx) && ffrac.collection{fridx,4}>40
            comidx=find(strcmp(fcom.collection(:,2),inter_reg(ri)));
            yy=fcom.collection{comidx,1}./4;
            xx=ffrac.collection{fridx,1};
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            %Notice, FRP still miss reg-color-group
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(inter_reg{ri},'large_area',true),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(inter_reg{ri},'large_area',true));
        end
    end
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);

    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    xlabel('Fraction of delay selective neuron')
    ylabel('F.R. center of mass (s)')
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    set(gca(),'XScale','linear');
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
%         exportgraphics(fh,sprintf('per_region_TCOM_FRAC_%d.pdf',opt.delay));
    end
end

%% com vs SMI
if opt.COM_sensemotor
%     load('OBM1Map.mat','OBM1map')
    fh=figure('Color','w','Position',[100,100,245,235]);
    hold on;
    coord=[];
    regs=[];
    inter_reg=intersect(fcom.collection(:,2),ureg);

    for ri=1:numel(inter_reg)
        if  hiermap.isKey(inter_reg{ri})
            comidx=find(strcmp(fcom.collection(:,2),inter_reg(ri)));
            xx=hiermap(inter_reg{ri});
            yy=fcom.collection{comidx,1}./4;
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(inter_reg{ri},'large_area',true),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(inter_reg{ri},'large_area',true));
        else
            warning('Missing PVSST map key');
            disp(inter_reg{ri})
        end
    end
    
    xlabel('Sensory Motor Index (A.U.)')
    ylabel('F.R. center of mass (s)')

    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    set(gca(),'XScale','log')
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');

    if opt.export
%         exportgraphics(fh,sprintf('per_region_COM_SMI_%d.pdf',opt.delay));
    end


end



