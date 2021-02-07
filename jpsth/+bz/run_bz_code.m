addpath('K:\refcode\buzcode\io')
addpath('K:\refcode\buzcode\utilities\')
addpath('K:\refcode\buzcode\analysis\spikes\correlation\')
addpath('K:\refcode\buzcode\visualization\')
basepath='K:\refcode\buzdata\20170301';
cd K:\refcode\buzcode\analysis\spikes\functionalConnectionIdentification
%% mono_res=bz_GetMonoSynapticallyConnected(basepath);

%load data
spikes = bz_GetSpikes('basepath',basepath);

%get shank clu
spikeIDs = [spikes.shankID(spikes.spindices(:,2))' spikes.cluID(spikes.spindices(:,2))' spikes.spindices(:,2)];

%call main script
mono_res = bz_MonoSynConvClick (spikeIDs,spikes.spindices(:,1),'plot',false);
