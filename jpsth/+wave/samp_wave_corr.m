% CONST
meta=ephys.util.load_meta();
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));

BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));

cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
grey_regs=grey_regs(cnt>100);


idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
load('map_cells.mat','map_cells');
%%
figure('Color','w')
for rii=1:numel(grey_regs)
    subplot(6,6,rii) 
    hold on;
    samp_dur_ep=nan(6,2);
    for epii=2:3
        sampep=map_cells{epii,1}(grey_regs{rii});
        durep=map_cells{epii,3}(grey_regs{rii});
        samp_dur_ep(epii,:)=[sampep(1),durep(1)];
    end
    plot(samp_dur_ep(:,1),samp_dur_ep(:,2),'-r')
    plot(samp_dur_ep(1,1),samp_dur_ep(1,2),'ko')
    tree=idmap.reg2tree(grey_regs{rii});
    title(tree([5,7]))
%     keyboard
end
