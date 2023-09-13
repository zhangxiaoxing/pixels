function [cover_per_sec,fh]=delay_vs_iti_per_sec(covered,trials_dict,opt)
arguments
    covered = []
    trials_dict = []
    opt.skip_save (1,1) logical = false
end
if isempty(covered)
    fstr=load(fullfile('binary','delay_iti_runlength_covered.mat'),'covered_tbl');
    covered=fstr.covered_tbl;
    clear fstr;
end
if isempty(trials_dict)
    load(fullfile('binary','trials_dict.mat'),'trials_dict');
end


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
    np_trl_sel=find(trials(:,5)==setdiff([4,8],samp) & trials(:,8)==delay & all(trials(:,9:10)>0,2));

    pref_delay_sec=sum(diff(trials(trl_sel,1:2),1,2)./sps-1);
    delayc=covered.delay{sessidx};
    delay_motif_sec=0;
    for tt=reshape(trl_sel,1,[])
        delay_motif_sec=delay_motif_sec+nnz(delayc(floor(trials(tt,1)./3):ceil(trials(tt,2)./3)))./10000;
    end
        
    np_delay_sec=sum(diff(trials(np_trl_sel,1:2),1,2)./sps-1);
    npdelayc=covered.npdelay{sessidx};
    npdelay_motif_sec=0;
    for tt=reshape(np_trl_sel,1,[])
        npdelay_motif_sec=npdelay_motif_sec+nnz(npdelayc(floor(trials(tt,1)./3):ceil(trials(tt,2)./3)))./10000;
    end

    trials(end+1,:)=trials(end,2)+14*sps;
    pref_succeed_iti_sec=sum((trials(trl_sel+1,1)-trials(trl_sel,2))./sps-4); % 1s test + 3s response
    itic=covered.iti{sessidx};
    iti_motif_sec=0;
    for tt=reshape(trl_sel,1,[])
        iti_motif_sec=iti_motif_sec+nnz(itic(floor(trials(tt,2)./3):ceil(trials(tt+1,1)./3)))./10000;
    end

    npiti_sec=sum((trials(np_trl_sel+1,1)-trials(np_trl_sel,2))./sps-4); % 1s test + 3s response
    npitic=covered.npiti{sessidx};
    npiti_motif_sec=0;
    for tt=reshape(np_trl_sel,1,[])
        npiti_motif_sec=npiti_motif_sec+nnz(npitic(floor(trials(tt,2)./3):ceil(trials(tt+1,1)./3)))./10000;
    end

    trials(end,:)=[];

    cover_per_sec=[cover_per_sec;...
        cell2table({sess,samp,delay,[delay_motif_sec,pref_delay_sec],[iti_motif_sec,pref_succeed_iti_sec],[surround_motif_sec,surround_sec],[npdelay_motif_sec,np_delay_sec],[npiti_motif_sec,npiti_sec]},...
        'VariableNames',{'session','sample','dur','delay','iti','surround','npdelay','npiti'})];

end

%% plot

dps=cover_per_sec.delay;
dprop=(dps(:,1)./dps(:,2));

npdps=cover_per_sec.npdelay;
npdprop=(npdps(:,1)./npdps(:,2));

ips=cover_per_sec.iti;
iprop=(ips(:,1)./ips(:,2));

npips=cover_per_sec.npiti;
npiprop=(npips(:,1)./npips(:,2));

sps=cover_per_sec.surround;
sprop=(sps(:,1)./sps(:,2)); % Due to 4 waves

mdm=median([dprop,npdprop,iprop,npiprop,sprop]);
ci=bootci(1000,@(x) median(x),[dprop,npdprop,iprop,npiprop,sprop]);
% mm=[mean(dprop),mean(npdprop),mean(iprop),mean(npiprop),mean(sprop),];
% sem=[std(dprop),std(npdprop),std(iprop),std(npiprop),std(sprop)]./sqrt([numel(dprop),numel(npdprop),numel(iprop),numel(npiprop),numel(sprop)]);

fh=figure('Position',[100,100,400,300]);
hold on
bh=bar(mdm.*100,'FaceColor','none','FaceColor','k');
errorbar(bh.XEndPoints,bh.YEndPoints,(ci(1,:)-mdm).*100,(ci(2,:)-mdm).*100,'k.')
set(gca,'XTick',1:5,'XTickLabelRotation',90,'XTickLabel',{'Delay','NPDelay','ITI','NPITI','Surround'});
ylabel('Motif duration / total duration (%)');

pnpdelay=signrank(dprop,npdprop);
piti=signrank(dprop,iprop);
piti_iti=signrank(iprop,npiprop);
psur=signrank(dprop,sprop);

title(sprintf('%s%.4f,','p-np',pnpdelay,'d-iti',piti,'iti-iti',piti_iti,'otask',psur));
savefig(fh,fullfile('binary','delay_vs_iti_per_sec.fig'));
appendfig('close',true,'tag','delay vs iti composite motif coverage time, delay_vs_iti_per_sec.m');
blame=vcs.blame();
save(fullfile('binary','motif_cover_per_sec.mat'),'cover_per_sec','blame');

end
