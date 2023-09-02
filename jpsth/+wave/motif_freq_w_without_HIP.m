function motif_freq_w_without_HIP()

load(fullfile('binary','motif_replay.mat'))

fhc=wave.replay.region_replay(chain_replay,'reg',"HIP");
ylim([0,1.5])
title('Chains, HIP')
% wave.replay.region_replay(chain_replay,'reg',"ORB")
% title('Chains, ORB')
fhl=wave.replay.region_replay(ring_replay,'reg',"HIP");
title('Loops, HIP')
ylim([0,1.5])
% wave.replay.region_replay(ring_replay,'reg',"ORB")
% title('Loops, ORB')
% ylim([0,1.5])

savefig(fhc,fullfile('binary','chain_freq_w_without_HIP.fig'))
savefig(fhl,fullfile('binary','loop_freq_w_without_HIP.fig'))