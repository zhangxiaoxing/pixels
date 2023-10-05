function motif_prefer_nonprefer()
load(fullfile('binary','motif_replay.mat'))
  % 12×1 cell array
  % 
  % 1  {'pref_delay_correct'   }
  % 2  {'pref_delay_error'     }
  % 3  {'nonpref_delay_correct'}
  % 4  {'pref_test'            }
  % 5  {'pref_succeed_ITI'     }
  % 6  {'pref_succeed_ITI_err' }
  % 7  {'nonpref_succeed_ITI'  }
  % 8  {'pref_precede_ITI'     }
  % 9  {'pref_precede_ITI_err' }
  % 10  {'nonpref_precede_ITI'  }
  % 11  {'before_session'       }
  % 12  {'after_session'        }
  % 13 {nonpref_delay_error}
% correct error chain
chain_corr_err=cell2struct({chain_raw.count+eps;...
    chain_raw.time;...
    chain_raw.condition; ...
    chain_raw.tag},{'count';'time';'condition';'tag'});

% correct error loop
ring_corr_err=cell2struct({loops_raw.count;...
    loops_raw.time;...
    loops_raw.condition; ...
    loops_raw.tag},{'count';'time';'condition';'tag'});

[~,jemat]=wave.replay.stats_replay_sess({chain_corr_err,ring_corr_err});
if false
cestr=cell2struct({jemat(:,1);jemat(:,2);jemat(:,3);jemat(:,13);jemat(:,8);jemat(:,9);jemat(:,5);jemat(:,6)},...
    {'preferred_delay_correct','preferred_delay_error','nonpreferred_delay_correct','nonpreferred_delay_error',...
    'ITI_before_correct_trials','ITI_before_error_trials','ITI_after_correct_trials','ITI_after_error_trials'});
fid=fopen(fullfile('binary','upload','F3H_Motif_spike_frequency_correct_error_trials.json'),'w');
fprintf(fid,jsonencode(cestr));
fclose(fid);
end


ratiomat=jemat(:,[2,13,9,6])./jemat(:,[1,3,8,5]); %npdelay/pdelay, npiti/piti
qtrs=prctile(ratiomat,[25,50,75]); % only median were used
bci=bootci(1000,@(x) median(x),ratiomat);

fce=figure('Position',[100,100,800,300]);
hold on
bh=bar([ones(1,4);qtrs(2,:)].','grouped');
bh(1).FaceColor='k';
bh(2).FaceColor='w';
errorbar(bh(2).XEndPoints,qtrs(2,:),bci(1,:)-qtrs(2,:),bci(2,:)-qtrs(2,:),'k.')
set(gca(),'XTick',1:4,'XTickLabel',{'Corr-err-delay','Corr-err-npdelay','ITI pre-corr-err','ITI-post-corr-err'})

pdelay=signrank(jemat(:,2),jemat(:,1));
pnpdelay=signrank(jemat(:,13),jemat(:,3));
ppre=signrank(jemat(:,9),jemat(:,8));
ppost=signrank(jemat(:,6),jemat(:,5));
title(sprintf('Correct/error,delay,npdelay,preITI,postITI,%.4f,%.4f,%.4f,%.4f',pdelay,pnpdelay,ppre,ppost))


ratiomat=jemat(:,[3,10,7])./jemat(:,[1,8,5]); %npdelay/pdelay, npiti/piti

if false
pnpstr=cell2struct({jemat(:,1);jemat(:,3);jemat(:,8);jemat(:,10);jemat(:,5);jemat(:,7)},...
    {'preferred_delay','nonpreferred_delay','ITI_before_preferred_trials','ITI_before_nonpreferred_trials',...
    'ITI_after_preferred_trials','ITI_after_nonpreferred_trials'});
fid=fopen(fullfile('binary','upload','SF7A_Motif_spike_frequency_preferred_nonpreferred_trials.json'),'w');
fprintf(fid,jsonencode(cestr));
fclose(fid);
end

qtrs=prctile(ratiomat,[25,50,75]);
bci=bootci(1000,@(x) median(x),ratiomat);

fpnp=figure('Position',[100,100,800,300]);
hold on
bh=bar([1,qtrs(2,1);1,qtrs(2,2);1,qtrs(2,3)],'grouped');
bh(1).FaceColor='k';
bh(2).FaceColor='w';
errorbar(bh(2).XEndPoints,qtrs(2,:),qtrs(1,:)-qtrs(2,:),qtrs(3,:)-qtrs(2,:),'k.')
set(gca(),'XTick',1:3,'XTickLabel',{'Prefer-npref-delay','ITI pre-pref-np','ITI-post-pref-np'})

pdelay=signrank(jemat(:,3),jemat(:,1));
ppre=signrank(jemat(:,10),jemat(:,8));
ppost=signrank(jemat(:,7),jemat(:,5));
title(sprintf('pref/nonpref,delay,preITI,postITI,%.4f,%.4f,%.4f',pdelay,ppre,ppost))

savefig(fce,fullfile("binary","motif_freq_correct_error.fig"));
savefig(fpnp,fullfile("binary","motif_freq_prefer_nonpref.fig"));

% 
% 
% [~,jemat]=wave.replay.stats_replay_sess({chain_corr_err,ring_corr_err});
% qtrs=prctile(jemat(:,[1,5,11,12]),[25,50,75]);
% fh=figure('Position',[100,100,800,300]);
% hold on
% bh=bar(qtrs(2,:),'grouped','FaceColor','k');
% 
% errorbar(bh.XEndPoints,qtrs(2,:),qtrs(1,:)-qtrs(2,:),qtrs(3,:)-qtrs(2,:),'k.')
% set(gca(),'XTick',1:3,'XTickLabel',{'Corr-err-delay','ITI pre-corr-err','ITI-post-corr-err'})
% 
% pdelay=signrank(jemat(:,2),jemat(:,1));
% ppre=signrank(jemat(:,9),jemat(:,8));
% ppost=signrank(jemat(:,6),jemat(:,5));
% title(sprintf('Correct/error,delay,preITI,postITI,%.4f,%.4f,%.4f',pdelay,ppre,ppost))




end
