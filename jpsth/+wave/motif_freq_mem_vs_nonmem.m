function [fh,cyy,cgg]=motif_freq_mem_vs_nonmem(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'chain','loop'})}= 'chain'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end

% TODO: Leanring nonmemory
if strcmp(opt.type,'chain')
    switch opt.criteria
        case 'WT'
            memstr=load(fullfile("binary","motif_replay.mat"),'chain_sums');
            nmstr=load(fullfile('binary','motif_replay_chain_nonmem.mat'),'chain_sums');
        case 'Learning'
            memstr=load(fullfile("binary","LN_motif_replay.mat"),'chain_sums');
            nmstr=load(fullfile('binary','LN_motif_replay_chain_nonmem.mat'),'chain_sums');            
        otherwise
            error("unfinished")
    end
    if false
        mvsnm.Memory=cell2struct(mat2cell(memstr.chain_sums([1 3 5 11 12],:),ones(5,1),size(memstr.chain_sums,2)),{'Preferred_delay','Nonpreferred_delay','ITI','Before_task','After_task'});
        mvsnm.Nonmemory=cell2struct(mat2cell(nmstr.chain_sums([1 1 2 4 5],:),ones(5,1),size(nmstr.chain_sums,2)),{'Preferred_delay','Nonpreferred_delay','ITI','Before_task','After_task'});
        fid=fopen(fullfile('binary','upload','F3F_Chains_frequency_memory_nonmemory.json'),'w');
        fprintf(fid,jsonencode(mvsnm));
        fclose(fid)
    end

    cyy=[memstr.chain_sums(1,:),nmstr.chain_sums(1,:),...
        memstr.chain_sums(3,:),nmstr.chain_sums(1,:),...
        memstr.chain_sums(5,:),nmstr.chain_sums(2,:),...
        memstr.chain_sums(11,:),nmstr.chain_sums(4,:),...
        memstr.chain_sums(12,:),nmstr.chain_sums(5,:)];
    ggn=[size(memstr.chain_sums,2),size(nmstr.chain_sums,2)];
else
    switch opt.criteria
        case 'WT'
            memstr=load(fullfile("binary","motif_replay.mat"),'loops_sums');
            nmstr=load(fullfile('binary','motif_replay_ring_nonmem.mat'),'loops_sums');
        case 'Learning'
            memstr=load(fullfile("binary","LN_motif_replay.mat"),'loops_sums');
            nmstr=load(fullfile('binary','LN_motif_replay_ring_nonmem.mat'),'loops_sums');
        otherwise
            error("Unfinished")
    end
    if false
        mvsnm.Memory=cell2struct(mat2cell(memstr.loops_sums([1 3 5 11 12],:),ones(5,1),size(memstr.loops_sums,2)),{'Preferred_delay','Nonpreferred_delay','ITI','Before_task','After_task'});
        mvsnm.Nonmemory=cell2struct(mat2cell(nmstr.loops_sums([1 1 2 4 5],:),ones(5,1),size(nmstr.loops_sums,2)),{'Preferred_delay','Nonpreferred_delay','ITI','Before_task','After_task'});
        fid=fopen(fullfile('binary','upload','F3G_Loops_frequency_memory_nonmemory.json'),'w');
        fprintf(fid,jsonencode(mvsnm));
        fclose(fid)
    end


    cyy=[memstr.loops_sums(1,:),nmstr.loops_sums(1,:),...
        memstr.loops_sums(3,:),nmstr.loops_sums(1,:),...
        memstr.loops_sums(5,:),nmstr.loops_sums(2,:),...
        memstr.loops_sums(11,:),nmstr.loops_sums(4,:),...
        memstr.loops_sums(12,:),nmstr.loops_sums(5,:)];
    ggn=[size(memstr.loops_sums,2),size(nmstr.loops_sums,2)];
end

cgg=[ones(ggn(1),1);2*ones(ggn(2),1);...
    3*ones(ggn(1),1);4*ones(ggn(2),1);...
    5*ones(ggn(1),1);6*ones(ggn(2),1);...
    7*ones(ggn(1),1);8*ones(ggn(2),1);...
    9*ones(ggn(1),1);10*ones(ggn(2),1)];

cmm=arrayfun(@(x) median(cyy(cgg==x & isfinite(cyy.'))),1:10);
cci=cell2mat(arrayfun(@(x) bootci(1000,@(x) median(x), cyy(cgg==x & isfinite(cyy.'))),1:10,'UniformOutput',false));

fh=figure();
hold on
bh=bar(cmm([1:2:10;2:2:10].'),'grouped','FaceColor','none','EdgeColor','k');
bh(1).FaceColor='k';
bh(2).FaceColor='w';
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,cci(1,1:2:10)-cmm(1:2:10),cci(2,1:2:10)-cmm(1:2:10),'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,cci(1,2:2:10)-cmm(2:2:10),cci(2,2:2:10)-cmm(2:2:10),'k.');
set(gca(),'XTick',1:5,'XTickLabel',{'Delay','NPDelay','ITI','Before','After'})

pp=[ranksum(cyy(cgg==1),cyy(cgg==2)),ranksum(cyy(cgg==3),cyy(cgg==4)),...
    ranksum(cyy(cgg==5),cyy(cgg==6)),ranksum(cyy(cgg==7),cyy(cgg==8)),...
    ranksum(cyy(cgg==9),cyy(cgg==10))];

subtitle(sprintf('%.4f,',pp));
switch opt.criteria
    case 'WT'
        title([opt.type,'s consis vs nonmem sample'])
        savefig(fh,fullfile('binary',['motif_freq_',opt.type,'_mem_vs_nonmem.fig']))
    case 'Learning'
        title(['LN_',opt.type,'s consis vs nonmem sample'])
        savefig(fh,fullfile('binary',['LN_motif_freq_',opt.type,'_mem_vs_nonmem.fig']))
    otherwise
        error("unfinished")
end


