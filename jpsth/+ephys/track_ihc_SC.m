%% find file
trkfile=["42_20191125_g0\42_20191125_g0_imec0_cleaned",...
    "42_20191126_g0\42_20191126_g0_imec0_cleaned",...
    "M42_20191118_g0\M42_20191118_g0_imec0_cleaned",...
    "M42_20191120_g0\M42_20191120_g0_imec0_cleaned",...
    "M42_20191121_g0\M42_20191121_g0_imec0_cleaned",...
    "M42_20191123_gc\M42_20191123_gc_imec0_cleaned"];

homedir=ephys.util.getHomedir('type','raw');

for dd=trkfile
    f=fullfile(homedir,dd,'FR_All_ 250.hdf5');
    if ~exist(f),continue;end
%     s=h5info(f);
%     disp({s.Datasets.Name});
    regtbl=readtable(fullfile(homedir,dd,'su_id2reg.csv'));
    disp(dd)
    disp('=======imec0========')
    disp(unique(regtbl.d7(regtbl.index<10000)))
    disp('=======imec1========')
    disp(unique(regtbl.d7(regtbl.index>=10000)))
    continue
end

%     42_20191126_g0\42_20191126_g0_imec0_cleaned
% =======imec0========
%     {'AI' }
%     {'AON'}
%     {'MO' }
%     {'ORB'}

%% plot
f=fullfile(homedir,'42_20191126_g0\42_20191126_g0_imec0_cleaned','FR_All_1000.hdf5');
metatbl=readtable(fullfile(homedir,'42_20191126_g0\42_20191126_g0_imec0_cleaned','cluster_info.tsv'),'FileType','text');
depmap=containers.Map(metatbl.id,metatbl.depth);

FR=h5read(f,'/FR_All');
trials=h5read(f,'/Trials');
suids=h5read(f,'/SU_id');
susel=suids<10000;
s1sel=trials(:,5)==4 & trials(:,8)==6 & all(trials(:,9:10),2);
s2sel=trials(:,5)==8 & trials(:,8)==6 & all(trials(:,9:10),2);

sudepth=[suids(susel),arrayfun(@(x) depmap(x),suids(susel))];

dep_sel_mat=nan(max(sudepth(:,2))./20,size(FR,3));
dense_mat=zeros(max(sudepth(:,2))./20,1);
for dd=20:20:max(sudepth(:,2))
    depsel=sudepth(:,2)==dd;
    dense_mat(dd/20)=nnz(depsel);
    if nnz(depsel)==0, continue;end
    frs1=squeeze(mean(FR(s1sel,depsel,:)));
    frs2=squeeze(mean(FR(s2sel,depsel,:)));
    if nnz(depsel)>1
        selidx=(mean(abs(frs1-frs2)./(frs1+frs2)));
    else
        selidx=(abs(frs1-frs2)./(frs1+frs2));
    end
    dep_sel_mat(dd/20,:)=selidx;
end


fh=figure('Color','w','Position',[100,100,150,240]);
hold on
colormap('parula')
imAlpha=ones(size(dep_sel_mat));
imAlpha(isnan(dep_sel_mat))=0;
im=imagesc(dep_sel_mat,[0,0.3]);
im.AlphaData=imAlpha;

arrayfun(@(x) xline(x,'--w'),[3,4,10,11]+0.5)
xlim([2.5,11.5])
xlabel('Time (s)')
ylabel('Distance to tip (mm)')

exportgraphics(fh,'track_ihc_selecIdx.pdf','ContentType','vector');

%% dense
densemat2=[dense_mat,dense_mat];
fh=figure('Color','w','Position',[100,100,150,240]);
hold on
colormap('parula')
imAlpha=ones(size(densemat2));
imAlpha(densemat2==0)=0;
im=imagesc(densemat2,[0,4]);
im.AlphaData=imAlpha;
xlim([0.5,2.5])
set(gca,'color',0*[1 1 1],'YDir','normal','YTick',0:50:150,'YTickLabel',0:3,'XTick',[]);
colorbar()

exportgraphics(fh,'track_ihc_density.pdf','ContentType','vector');



%% reward

psel=trials(:,5)~=trials(:,6) & trials(:,8)==6 & all(trials(:,9:10),2);
npsel=trials(:,5)==trials(:,6) & trials(:,8)==6 & all(trials(:,9:10),2);

for dd=20:20:max(sudepth(:,2))
    depsel=sudepth(:,2)==dd;
    if nnz(depsel)==0, continue;end
    frs1=squeeze(mean(FR(psel,depsel,:)));
    frs2=squeeze(mean(FR(npsel,depsel,:)));
    if nnz(depsel)>1
        selidx=(mean(abs(frs1-frs2)./(frs1+frs2)));
    else
        selidx=(abs(frs1-frs2)./(frs1+frs2));
    end
    dep_sel_mat(dd/20,:)=selidx;
end

fh=figure('Color','w','Position',[100,100,150,240]);
hold on
colormap('parula')
imAlpha=ones(size(dep_sel_mat));
imAlpha(isnan(dep_sel_mat))=0;
im=imagesc(dep_sel_mat,[0.,0.3]);
im.AlphaData=imAlpha;
set(gca,'color',0*[1 1 1],'YDir','normal','YTick',0:25:150,'YTickLabel',0:0.5:3,'XTick',[4.5,9.5,14.5],'XTickLabel',[0 5 10]);
arrayfun(@(x) xline(x,'--w'),[3,4,10,11]+0.5)
xlim([9.5,14.5])
xlabel('Time (s)')
ylabel('Distance to tip (mm)')

exportgraphics(fh,'track_ihc_rewardIdx.pdf','ContentType','vector');

fh=figure('Color','w','Position',[100,100,150,240]);
hold on
colormap('parula')
imAlpha=ones(size(dep_sel_mat));
imAlpha(isnan(dep_sel_mat))=0;
im=imagesc(dep_sel_mat,[0.,0.3]);
im.AlphaData=imAlpha;
set(gca,'color',0*[1 1 1],'YDir','normal','YTick',0:25:150,'YTickLabel',0:0.5:3,'XTick',[4.5,9.5,14.5],'XTickLabel',[0 5 10]);
arrayfun(@(x) xline(x,'--w'),[3,4,10,11]+0.5)
xlim([9.5,14.5])
xlabel('Time (s)')
ylabel('Distance to tip (mm)')
colorbar()

exportgraphics(fh,'track_ihc_rewardIdx_colorbar.pdf','ContentType','vector');
