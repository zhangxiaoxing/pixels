meta=ephys.util.load_meta();
waveid=ephys.get_wave_id(meta.sess,meta.allcid);

% mean(strcmp(meta.reg_tree(4,waveid>4),'OLF'))
% mean(strcmp(meta.reg_tree(4,waveid>0 & waveid<5),'OLF'))

load('OBM1map.mat','OBM1map');
regremat=[];

figure('Color','w')
hold on
for regs=OBM1map.keys
    mapv=OBM1map(regs{1});
    both=nnz(strcmp(meta.reg_tree(5,:),regs{1}).' & waveid>4);
    either=nnz(strcmp(meta.reg_tree(5,:),regs{1}).' & waveid>0 & waveid<5);
    bothrate=both./(both+either).*100;
    plot(mapv,bothrate,'r.')
    text(mapv,bothrate,regs{1},'HorizontalAlignment','center','VerticalAlignment','top')
    regremat=[regremat;mapv,bothrate];
end
[rr,pp]=corr(regremat(:,1),regremat(:,2));
regremat(:,3)=1;
regres=regremat(:,[1,3])\regremat(:,2);
plot(xlim(),xlim().*regres(1)+regres(2),'--k');
text(max(xlim()),max(ylim()),sprintf('r=%0.3f,p=%0.3f',rr,pp),'HorizontalAlignment','right','VerticalAlignment','top');

xlabel('Olfactory-motor connectivity index')
ylabel('Both-duration selective / total selective (%)')