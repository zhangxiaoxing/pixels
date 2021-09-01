keyboard();
sess=18;
fs=sprintf('s%d',sess);
com_map=wave.get_com_map('onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',true,'rnd_half',true);

samp_key_S1=intersect(cell2mat(com_map.(fs).s1a.keys),cell2mat(com_map.(fs).s1b.keys));
samp_key_S2=intersect(cell2mat(com_map.(fs).s2a.keys),cell2mat(com_map.(fs).s2b.keys));
COMS1=cell2mat(values(com_map.(fs).s1a,num2cell(samp_key_S1)));
COMS2=cell2mat(values(com_map.(fs).s2a,num2cell(samp_key_S2)));

r=corr([COMS1,COMS2].',...
    [cell2mat(values(com_map.(fs).s1b,num2cell(samp_key_S1))),...
    cell2mat(values(com_map.(fs).s2b,num2cell(samp_key_S2)))].');
sortmat=[ones(size(COMS1)),2*ones(size(COMS2));...
    double(samp_key_S1),double(samp_key_S2);...
    COMS1,COMS2].';

sortmat=sortrows(sortmat,3);

immata=[];
for ri=1:size(sortmat,1)
    one_heat=com_map.(fs).(sprintf('s%daheat',sortmat(ri,1)))(sortmat(ri,2));
    immata=[immata;one_heat];
end
immatb=[];
com_b=[];
for ri=1:size(sortmat,1)
    one_heat=com_map.(fs).(sprintf('s%dbheat',sortmat(ri,1)))(sortmat(ri,2));
    immatb=[immatb;one_heat];
    com_b=[com_b;com_map.(fs).(sprintf('s%db',sortmat(ri,1)))(sortmat(ri,2))];
end

fh=figure('Color','w','Position',[32,32,300,350]);
plotOne(1,immata,sortmat(:,3));
plotOne(2,immatb,com_b);
shufidx=randsample(size(immatb,1),size(immatb,1));
plotOne(3,immatb(shufidx,:),com_b(shufidx));
sgtitle(num2str(r));
exportgraphics(fh,'wave_half_half.pdf')

function plotOne(subidx,imdata,comdata)
subplot(3,1,subidx);
hold on
colormap('jet');
gk = fspecial('gaussian', [3 3], 1);
imagesc(conv2(imdata,gk,'same'),[-0.8 0.8])
scatter(comdata,1:numel(comdata),2,'o','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
colorbar();
ylim([0.5,size(imdata,1)+0.5])
xlim([0.5,24.5])
end