nestedloopfstr=load(fullfile('binary','delay_iti_runlength_covered.mat'),'run_length');
load(fullfile('binary','motif_replay.mat'),'chain_replay','ring_replay');
chainrunlength=[];
for cii=1:size(chain_replay,1)
    trl_align=chain_replay.trl_align{cii};
    pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
    chainrunlength=[chainrunlength;arrayfun(@(x) diff(chain_replay.ts{cii}(x,[1 end]),1,2), find(pref_delay))./30];
end
looprunlength=[];
for rii=1:size(ring_replay,1)
    trl_align=ring_replay.trl_align{rii};
    pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
    looprunlength=[looprunlength;cellfun(@(x) x(end)-x(1), ring_replay.ts{rii}(pref_delay))./30];
end

chainhist=histcounts(chainrunlength,[15:5:45],'Normalization','pdf');
chainf=fit((17.5:5:42.5).',chainhist.','power1','MaxIter',4000,Lower=[0,-inf],Upper=[inf,0],Weights=1./chainhist)

loophist=histcounts(looprunlength,[25:10:85],'Normalization','pdf');
loopf=fit((27.5:10:82.5).',loophist.','power1','MaxIter',4000,Lower=[0,-inf],Upper=[inf,0],Weights=1./loophist)

nestedloopspdf=histcounts(loopfstr.run_length.delay(:,2),[10:20:280],'Normalization','pdf');
nestf=fit((20:20:270).',nestedloopspdf.','power1','MaxIter',4000,Lower=[0,-inf],Upper=[inf,0],Weights=1./nestedloopspdf)

coef=[coeffvalues(chainf),coeffvalues(loopf),coeffvalues(nestf)];

figure()
set(gca(),'XScale','log','YScale','log')
hold on
plot(17.5:5:42.5,chainhist,'b:')
plot(chainf,'b-')
plot((27.5:10:82.5),loophist,'r:')
plot((27.5:10:82.5),feval(loopf,(27.5:10:82.5)),'r-')
plot(20:20:270,nestedloopspdf,'k:')
plot(nestf,'k-')
xlim([5,5e2])
legend('off')
xlabel('Time (ms)')
ylabel('PDF')
title(sprintf('power law v time const c%.3fl%.3fn%.3f',coef([2 4 6])))

% modelfun=@(b,x) b(1)*(x.^b(2)+b(3));
% xx=((12.5:5:50)).';
% tbl=table(xx,chainhist().','VariableNames',{'time','prob'});
% chainmdl=fitnlm(tbl,modelfun,[0.5,0.5,0],'Options',statset('MaxIter',1000));
% xx=(25:10:90).';
% tbl=table(xx,loophist.','VariableNames',{'time','prob'});
% loopmdl=fitnlm(tbl,modelfun,[0.5,0.5,0],'Options',statset('MaxIter',1000));
% xx=(30:20:320).';
% tbl=table(xx,nestedloopspdf.','VariableNames',{'time','prob'});
% nestmdl=fitnlm(tbl,modelfun,[0.5,0.5,0],'Options',statset('MaxIter',1000));





% plot(nestf,'k-')

% plot(loopf,'r-')

% plot(nestmdl.Variables.time,nestmdl.Fitted,'-k','LineWidth',1);
% plot(chainmdl.Variables.time,chainmdl.Fitted,'-b','LineWidth',1)
% plot(loopmdl.Variables.time,loopmdl.Fitted,'-r','LineWidth',1)