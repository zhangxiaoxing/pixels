[~,samepair]=wave.com_diff('level',5,'diff',false,'to_plot',false,'per_sec_stats',false);
[~,diffpair]=wave.com_diff('level',5,'diff',true,'to_plot',false,'per_sec_stats',false);

same_dir=samepair(samepair(:,2)>samepair(:,1),:)./1000;
diff_dir=diffpair(diffpair(:,2)>diffpair(:,1),:)./1000;
[~,sameidx]=sort(mean(same_dir.'));
[~,diffidx]=sort(mean(diff_dir.'));

fh=figure('Color','w','Position',[100,100,235,235]);
subplot(1,2,1)
hold on;
plot(same_dir(sameidx,:).',ones(2,1)*(1:size(same_dir)),'-k');
ylim([0,size(same_dir,1)])
xlim([0,6]);
xlabel('Delay Time (s)')
title('Within region')
set(gca(),'FontSize',10)
subplot(1,2,2)
hold on;
plot(diff_dir(diffidx,:).',ones(2,1)*(1:size(diff_dir)),'-r');
ylim([0,size(diff_dir,1)]);
xlim([0,6]);
xlabel('Delay Time (s)')
title('Cross region')
set(gca(),'FontSize',10)
% exportgraphics(fh,'COM_pre_post.pdf','ContentType','vector');
