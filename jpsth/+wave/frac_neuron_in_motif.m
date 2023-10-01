function frac_neuron_in_motif()
% arguments
%     opt.shuf (1,1) boolean = false
%     opt.shufidx (1,2) double =[1,100]
% end


load(fullfile('binary','motif_replay.mat'),'chain_replay','ring_replay');
load(fullfile('binary','wrs_mux_meta.mat'),'wrs_mux_meta');
memory_su_cnt=nnz(wrs_mux_meta.wave_id>0 & wrs_mux_meta.wave_id<7);
allsu=cell2mat(...
    [arrayfun(@(x) [repmat(chain_replay.session(x),numel(chain_replay.meta{x,2}),1),chain_replay.meta{x,2}.'],(1:size(chain_replay,1)).','UniformOutput',false);...
    arrayfun(@(x) int32([repmat(ring_replay.session(x),numel(ring_replay.meta{x,2}),1),ring_replay.meta{x,2}.']),(1:size(ring_replay,1)).','UniformOutput',false)]);
motif_su_cnt=size(unique(allsu,'rows'),1);

[bhat,bci]=binofit(motif_su_cnt,memory_su_cnt);
shufcnt=zeros(100,1);
for shufidx=1:100
    load(fullfile("binary","motif_replay_shuf"+shufidx+".mat"),'chain_replay','ring_replay');
    shufsu=cell2mat(...
        [arrayfun(@(x) [repmat(chain_replay.session(x),numel(chain_replay.meta{x,2}),1),chain_replay.meta{x,2}.'],(1:size(chain_replay,1)).','UniformOutput',false);...
        arrayfun(@(x) [repmat(ring_replay.session(x),numel(ring_replay.meta{x,2}),1),ring_replay.meta{x,2}.'],(1:size(ring_replay,1)).','UniformOutput',false)]);
    shufcnt(shufidx)=size(unique(shufsu,'rows'),1);
end
[shat,sci]=binofit(sum(shufcnt),numel(shufcnt)*memory_su_cnt);

if false
    % usu=array2table(unique(allsu,'rows'),...
    %     'VariableNames',{'Session','ID'});
    % load(fullfile('binary','su_meta.mat'),'su_meta');
    % wmsu=array2table([su_meta.sess(wrs_mux_meta.wave_id>0 & wrs_mux_meta.wave_id<7),su_meta.allcid(wrs_mux_meta.wave_id>0 & wrs_mux_meta.wave_id<7)],...
    %     'VariableNames',{'Session','ID'});
    % motifsu=cell2struct({wmsu;usu},{'Memory_neuron','Motif_memory_neuron'});
    % fid=fopen(fullfile('binary','upload','F2N_motif_neuron_over_WM_neuron.json'),'w');
    % fprintf(fid,jsonencode(motifsu));
    % fclose(fid);
    motifsu.observed_motif_neuron_count=motif_su_cnt;
    motifsu.memory_neuron_count=memory_su_cnt;
    motifsu.shuffled_motif_neuron_count=shufcnt;

    fid=fopen(fullfile('binary','upload','F2N_motif_neuron_over_WM_neuron.json'),'w');
    fprintf(fid,jsonencode(motifsu));
    fclose(fid);

end




fh=figure();
hold on
bh=bar([bhat,shat].*eye(2),'stacked','FaceColor','w');
bh(1).FaceColor='k';
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,[bci(1),sci(1)]-[bhat,shat],[bci(2),sci(2)]-[bhat,shat],'k.','CapSize',15)
ylim([0,0.2]);
set(gca,'YTick',0:0.1:0.2,'YTickLabel',0:10:20,'XTick',[])
ylabel('WM neurons associated with motifs(%)')
% appendfig('tag','motif neuron / wm neuron, chain_plots.m')
legend(bh,{'Observed','Shuffled'})
savefig(fh,fullfile('binary','frac_neuron_in_motif.fig'))


[~,~,p]=crosstab([ones(1,memory_su_cnt),2*ones(1,100*memory_su_cnt)],[(1:memory_su_cnt)>motif_su_cnt,(1:(100*memory_su_cnt))>sum(shufcnt)]);