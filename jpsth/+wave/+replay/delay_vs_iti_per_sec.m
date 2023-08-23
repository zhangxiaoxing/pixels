function [cover_per_sec,fh]=delay_vs_iti_per_sec(covered,trials_dict)
sps=30000;
% per session
cover_per_sec=[];
for sessidx=1:size(covered,1)
    sess=covered.session(sessidx);
    switch covered.wave(sessidx)
        case "s1d3"
            samp=4;delay=3;
        case "s1d6"
            samp=4;delay=6;
        case "s2d3"
            samp=8;delay=3;
        case "s2d6"
            samp=8;delay=6;
    end

    disp(sess)
    session_tick=wave.replay.sessid2length(sess);
    trials=cell2mat(trials_dict(sess));
    surround_sec=(trials(1,1)./sps-60)+((session_tick-trials(end,2))./sps-60-1);
    surround_motif_sec=nnz(covered.out_task{sessidx})./10000;

    %  corresponding network in pre task, post task

    % per preferred trial
    trl_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
    pref_delay_sec=sum(diff(trials(trl_sel,1:2),1,2)./sps-1);
    delayc=covered.delay{sessidx};
    delay_motif_sec=0;
    for tt=reshape(trl_sel,1,[])
        delay_motif_sec=delay_motif_sec+nnz(delayc(floor(trials(tt,1)./3):ceil(trials(tt,2)./3)))./10000;
    end
    trials(end+1,:)=trials(end,2)+14*sps;
    pref_succeed_iti_sec=sum((trials(trl_sel+1,1)-trials(trl_sel,2))./sps-4); % 1s test + 3s response
    itic=covered.iti{sessidx};
    iti_motif_sec=0;
    for tt=reshape(trl_sel,1,[])
        iti_motif_sec=iti_motif_sec+nnz(itic(floor(trials(tt,2)./3):ceil(trials(tt+1,1)./3)))./10000;
    end
    trials(end,:)=[];

    cover_per_sec=[cover_per_sec;...
        cell2table({sess,samp,delay,[delay_motif_sec,pref_delay_sec],[iti_motif_sec,pref_succeed_iti_sec],[surround_motif_sec,surround_sec]},...
        'VariableNames',{'session','sample','dur','delay','iti','surround'})];

end

%% plot

dps=cell2mat(cover_per_sec.delay);
dprop=(dps(:,1)./dps(:,2));

ips=cell2mat(cover_per_sec.iti);
iprop=(ips(:,1)./ips(:,2));

sps=cell2mat(cover_per_sec.surround);
sprop=(sps(:,1)./sps(:,2))./4; % Due to 4 waves

mm=[mean(dprop),mean(iprop),mean(sprop)];
sem=[std(dprop),std(iprop),std(sprop)]./sqrt([numel(dprop),numel(iprop),numel(sprop)]);


fh=figure('Position',[100,100,400,300]);
hold on
bh=bar(mm.*100,'FaceColor','none');
errorbar(bh.XEndPoints,bh.YEndPoints,sem.*100,'k.')
set(gca,'XTick',1:3,'XTickLabelRotation',90,'XTickLabel',{'Delay','ITI','Surround'});
ylabel('Motif duration / total duration (%)');

[~,piti]=ttest(dprop,iprop);
[~,psur]=ttest2(dprop,sprop);

title(sprintf('%.4f,',piti,psur));

% keyboard()
end
