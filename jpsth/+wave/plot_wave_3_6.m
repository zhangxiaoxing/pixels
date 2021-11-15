keyboard();
for sess=102
%     while true
        fs=sprintf('s%d',sess);
        com_map_3=wave.get_com_map('onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',true,'rnd_half',false,'delay',3);
        com_map_6=wave.get_com_map('onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',true,'rnd_half',false,'delay',6,'early3in6',true);
        if ~isfield(com_map_6,fs) || ~isfield(com_map_3,fs)
            continue
        end

        samp_key_S1=intersect(cell2mat(com_map_3.(fs).s1.keys),cell2mat(com_map_6.(fs).s1.keys));
        samp_key_S2=intersect(cell2mat(com_map_3.(fs).s2.keys),cell2mat(com_map_6.(fs).s2.keys));
        COMS1=cell2mat(values(com_map_3.(fs).s1,num2cell(samp_key_S1)));
        COMS2=cell2mat(values(com_map_3.(fs).s2,num2cell(samp_key_S2)));
        if numel(COMS1)+numel(COMS2)<50
            continue;
        end
        r=corr([COMS1,COMS2].',...
            [cell2mat(values(com_map_6.(fs).s1,num2cell(samp_key_S1))),...
            cell2mat(values(com_map_6.(fs).s2,num2cell(samp_key_S2)))].');

%         if r<0.8
%             continue
%         end
%         keyboard()
        sortmat=[ones(size(COMS1)),2*ones(size(COMS2));...
            double(samp_key_S1),double(samp_key_S2);...
            COMS1,COMS2].';

        sortmat=sortrows(sortmat,3);

        immata=[];
        for ri=1:size(sortmat,1)
            one_heat=com_map_3.(fs).(sprintf('s%dcurve',sortmat(ri,1)))(sortmat(ri,2));
            immata=[immata;one_heat];
        end
        immatb=[];
        com_b=[];
        for ri=1:size(sortmat,1)
            one_heat=com_map_6.(fs).(sprintf('s%dcurve',sortmat(ri,1)))(sortmat(ri,2));
            immatb=[immatb;one_heat];
            com_b=[com_b;com_map_6.(fs).(sprintf('s%d',sortmat(ri,1)))(sortmat(ri,2))];
        end
        immata=immata./max(immata,[],2);
        immatb=immatb./max(immatb,[],2);

        fh=figure('Color','w','Position',[32,32,1080,140]);
        plotOne(1,immata,sortmat(:,3));
        plotOne(2,immatb,com_b);
        shufidx=randsample(size(immatb,1),size(immatb,1));
        plotOne(3,immatb(shufidx,:),com_b(shufidx));
        sgtitle(num2str([sess,r]));
%         exportgraphics(fh,sprintf('wave_half_half_%d.png',sess));
        keyboard()        
        exportgraphics(fh,sprintf('wave_3s_6s_%d.pdf',sess));
%         close all
%         waitfor(fh)
        
%     end
end

function plotOne(subidx,imdata,comdata)
subplot(1,3,subidx);
hold on
colormap('turbo');
gk = fspecial('gaussian', [3 3], 1);
imagesc(conv2(imdata,gk,'same'),[-1 1])
scatter(comdata,1:numel(comdata),2,'o','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
if size(imdata,2)>20
    set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
    xline(12.5,'--w','LineWidth',1);
else
    set(gca(),'XTick',[0.5,12.5],'XTickLabel',[0,3]);
end
colorbar();
ylim([0.5,size(imdata,1)+0.5])
xlim([0.5,size(imdata,2)+0.5])
end