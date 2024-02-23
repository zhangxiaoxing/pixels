function per_sess_bz_corr(sess_path,negccg)
arguments
	sess_path
	negccg (1,1) logical = false
end
disp(['processing ',sess_path])
ephys.util.dependency("buz",true,"ft",false)
addpath(fullfile("..","..","npy-matlab","npy-matlab"))
spkarr=readNPY(fullfile(sess_path,"sess_spk_t_id.npy"));
all_probe_ccg=[];
uprobe=unique(spkarr(:,1));
for ii=reshape(uprobe,1,[])
	prbsel=spkarr(:,1)==ii;
	mono=bz.sortSpikeIDz(spkarr(prbsel,2),spkarr(prbsel,3),'negccg',negccg,'wallclocktime',true,'binSize',0.0004); % adapted from English, Buzsaki, 2017
	mono.probe_id=ii;
	mono.sess_path=sess_path;
	all_probe_ccg=[all_probe_ccg;mono]
end

blame=vcs.blame();
if negccg
	fname='bz_corr_negccg.mat';
else
	fname='bz_corr.mat';
end	
save(fullfile(sess_path,fname),'all_probe_ccg','blame')
