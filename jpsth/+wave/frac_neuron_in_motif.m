function frac_neuron_in_motif()

load(fullfile('binary','motif_replay.mat'),'chain_replay','ring_replay');
load(fullfile('binary','wrs_mux_meta.mat'),'wrs_mux_meta');
if false
    allsu=cell2mat(...
        [arrayfun(@(x) [repmat(chain_replay.session(x),numel(chain_replay.meta{x,2}),1),chain_replay.meta{x,2}.'],(1:size(chain_replay,1)).','UniformOutput',false);...
        arrayfun(@(x) int32([repmat(ring_replay.session(x),numel(ring_replay.meta{x,2}),1),ring_replay.meta{x,2}.']),(1:size(ring_replay,1)).','UniformOutput',false)]);
    usu=array2table(unique(allsu,'rows'),...
        'VariableNames',{'Session','ID'});
    load(fullfile('binary','su_meta.mat'),'su_meta');
    wmsu=array2table([su_meta.sess(wrs_mux_meta.wave_id>0 & wrs_mux_meta.wave_id<7),su_meta.allcid(wrs_mux_meta.wave_id>0 & wrs_mux_meta.wave_id<7)],...
        'VariableNames',{'Session','ID'});
    motifsu=cell2struct({wmsu;usu},{'Memory_neuron','Motif_memory_neuron'});
    fid=fopen(fullfile('binary','upload','F2N_motif_neuron_over_WM_neuron.json'),'w');
    fprintf(fid,jsonencode(motifsu));
    fclose(fid);
end

[bhat,bci]=binofit(numel(unique([...
cell2mat(arrayfun(@(x) double(chain_replay.session(x)).*100000+double(chain_replay.meta{x,2}),1:size(chain_replay,1),'UniformOutput',false)),...
cell2mat(arrayfun(@(x) double(ring_replay.session(x)).*100000+double(ring_replay.meta{x,2}),1:size(ring_replay,1),'UniformOutput',false))...
])),nnz(wrs_mux_meta.wave_id>0 & wrs_mux_meta.wave_id<7));

fh=figure();
hold on
bar(bhat,'FaceColor','w')
errorbar(1,bhat,bci(1)-bhat,bci(2)-bhat,'k.')
ylim([0,0.2]);
set(gca,'YTick',0:0.1:0.2,'YTickLabel',0:10:20,'XTick',[])
% appendfig('tag','motif neuron / wm neuron, chain_plots.m')
savefig(fullfile('binary','frac_neuron_in_motif.fig'))