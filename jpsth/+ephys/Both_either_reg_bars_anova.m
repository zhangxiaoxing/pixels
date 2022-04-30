% pending remove? 
keyboard()
load('cross_ep_anovameta.mat','cross_ep_anovameta'); % K:\code\jpsth\+ephys\selectivity_anova.m
fn=fieldnames(cross_ep_anovameta);
meta=ephys.util.load_meta();
pdffile='reg_corr_anova.pdf';
for fii=1:numel(fn)
    anovameta=cross_ep_anovameta.(fn{fii});

    % waveid=ephys.get_wave_id(meta.sess,meta.allcid);

    mixed_sample_sel=any(anovameta.anovap(:,[1 3])<0.05,2);
    mixed_dur_sel=any(anovameta.anovap(:,2:3)<0.05,2);
    exclu_sample_sel=mixed_sample_sel & ~mixed_dur_sel;
    exclu_dur_sel=mixed_dur_sel & ~mixed_sample_sel;
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    ureg=unique(meta.reg_tree(5,strcmp(meta.reg_tree(1,:),'CH') | strcmp(meta.reg_tree(1,:),'BS')));
    sums=[];
    for reg=reshape(ureg,1,[])
        if isempty(reg{1}),continue;end
        regsel=strcmp(meta.reg_tree(5,:),reg).';
        cnt=nnz(regsel);
        if cnt>100
            bothcnt=nnz(regsel & mixed_sample_sel & mixed_dur_sel);
            exclu_sample_cnt=nnz(regsel & exclu_sample_sel);
            exclu_dur_cnt=nnz(regsel & exclu_dur_sel);
            mixed_sample_cnt=nnz(regsel & mixed_sample_sel);
            mixed_dur_cnt=nnz(regsel & mixed_dur_sel);

            grp=idmap.reg2tree(reg{1});
            sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,bothcnt,exclu_sample_cnt,exclu_dur_cnt,mixed_sample_cnt,mixed_dur_cnt];
            %===========================================================3====4==========5==============6==============7==================8=======
        end
    end

    sums(:,9:23)=0;

    for ii=1:size(sums,1)
        [phatb,  pcib]=binofit(sums(ii,4),sums(ii,3));
        [phatsx,pcisx]=binofit(sums(ii,5),sums(ii,3));
        [phatdx,pcidx]=binofit(sums(ii,6),sums(ii,3));
        [phatsm,pcism]=binofit(sums(ii,7),sums(ii,3));
        [phatdm,pcidm]=binofit(sums(ii,8),sums(ii,3));
        sums(ii,9:23)=[phatb,pcib,phatsx,pcisx,phatdx,pcidx,phatsm,pcism,phatdm,pcidm];
        %===============9===========12===========15============18===========21
    end


    fh=scatterOne(sums, 18,12,'Mixed sample selective proportion','Exclusive sample selective proportion',diag=true,ttl=fn{fii},pdffile=pdffile);
    fh=scatterOne(sums, 21,15,'Mixed duration selective proportion','Exclusive duration selective proportion',diag=true,ttl=fn{fii},pdffile=pdffile);
    fh=scatterOne(sums, 9,12,'Both sample selective proportion','Exclusive sample selective proportion',diag=true,ttl=fn{fii},pdffile=pdffile);
    fh=scatterOne(sums, 9,15,'Both duration selective proportion','Exclusive duration selective proportion',diag=true,ttl=fn{fii},pdffile=pdffile);
    fh=scatterOne(sums, 9,18,'Both sample selective proportion','Mixed sample selective proportion',diag=true,ttl=fn{fii},pdffile=pdffile);
    fh=scatterOne(sums, 9,21,'Both duration selective proportion','Mixed duration selective proportion',diag=true,ttl=fn{fii},pdffile=pdffile);
    fh=scatterOne(sums, 18,21,'Mixed sample selective proportion','Mixed duration selective proportion',ttl=fn{fii},pdffile=pdffile);
    fh=scatterOne(sums, 12,15,'Exclusive sample selective proportion','Exclusive duration selective proportion',ttl=fn{fii},pdffile=pdffile);

%     fh=barOne(sums,18,12,'Mixed sample','Exclusive sample');
end
% scatter >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function fh=scatterOne(sums,xcol,ycol,xlbl,ylbl,opt)
arguments
    sums
    xcol (1,1) double {mustBeInteger,mustBePositive}
    ycol (1,1) double {mustBeInteger,mustBePositive}
    xlbl (1,:) char
    ylbl (1,:) char
    opt.xlim = []
    opt.ylim = []
    opt.ttl = []
    opt.diag (1,1) logical = false
    opt.pdffile = []
end
persistent idmap flatten
if isempty(idmap) || isempty(flatten)
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    flatten=@(y) cellfun(@(x) x,y);
end

fh=figure('Color','w');
hold on
scatter(sums(:,xcol),sums(:,ycol))
text(sums(:,xcol),sums(:,ycol),flatten(idmap.ccfid2reg.values(num2cell(sums(:,2)))),'HorizontalAlignment','center','VerticalAlignment','top')
[r,p]=corr(sums(:,xcol),sums(:,ycol));

set(gca(),'XScale','linear','YScale','linear')

if ~isempty(opt.xlim), xlim(opt.xlim);end
if ~isempty(opt.ylim), ylim(opt.ylim);end

xlabel(xlbl)
ylabel(ylbl)
text(min(xlim()),max(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom','Color','r')
if opt.diag
    span=max([xlim(),ylim()]);
    plot([0,span],[0,span],'--k')
end
if ~isempty(opt.ttl),title(opt.ttl);end
if ~isempty(opt.pdffile)
    exportgraphics(fh,opt.pdffile,'ContentType','vector','Append',true);
end
end
%v==================================================
function fh=barOne(sums,cola,colb,lgda,lgdb,opt)
arguments
    sums
    cola (1,1) double {mustBeInteger,mustBePositive}
    colb (1,1) double {mustBeInteger,mustBePositive}
    lgda (1,:) char
    lgdb (1,:) char
    opt.ylim = []
end
persistent idmap flatten
if isempty(idmap) || isempty(flatten)
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    flatten=@(y) cellfun(@(x) x,y);
end

bardata=sortrows(sums,cola,'descend');

fh=figure('Color','w','Position',[32,32,1024,360]);
hold on
bh=bar(bardata(:,[cola,colb]),'grouped');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,cola+(0:1)),1,2),diff(bardata(:,cola+[0,2]),1,2),'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,colb+(0:1)),1,2),diff(bardata(:,colb+[0,2]),1,2),'k.');
bh(1).FaceColor='r';
bh(2).FaceColor='b';
set(gca(),'YScale','Linear')
if ~isempty(opt.ylim),ylim(opt.ylim);end
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2)))),'XTickLabelRotation',90)
% exportgraphics(fh,'Both_either_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))))),'UniformOutput',false)
ylabel('Proportion of selective neuron')
set(gca(),'YTick',0:0.1:max(opt.ylim),'YTickLabel',0:10:max(opt.ylim.*100))
legend(bh,{lgda,lgdb})
end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



