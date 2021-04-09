function mem_nonmem_dist(rings,opt)
% 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched
arguments
    rings (:,3) cell
    opt.subsel (1,1) logical = false
end

homedir = fullfile('K:','code','per_sec');
reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
cid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id');
allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
mem_type=h5read(fullfile(homedir,'transient_6.hdf5'),'/mem_type');

trial_sel=all(trial_counts>=20,1);
% wf_sel=(wf_good>0)';
reg_sel=strcmp(reg_tree(1,:),'BS') | strcmp(reg_tree(1,:),'CH');
subsel=trial_sel & reg_sel;
if opt.subsel
    cid=int32(cid(subsel));
    allpath=allpath(subsel);
    mem_type=mem_type(subsel)';
end
sess=arrayfun(@(x) ephys.path2sessid(allpath{x}),1:numel(allpath));
ucid=sess'.*100000+int32(cid);

sess_tagged=arrayfun(@(x) {x*100000+rings{x,1},...
    x*100000+rings{x,2},...
    x*100000+rings{x,3}},...
    1:size(rings,1),'UniformOutput',false);
[SUC,SU]=groupcounts([reshape(cell2mat(cellfun(@(x) x{1},sess_tagged,'UniformOutput',false)'),[],1);...
    reshape(cell2mat(cellfun(@(x) x{2},sess_tagged,'UniformOutput',false)'),[],1);...
    reshape(cell2mat(cellfun(@(x) x{3},sess_tagged,'UniformOutput',false)'),[],1)]);

sust_ucid=ucid(ismember(mem_type,[1,3]));
trans_ucid=ucid(ismember(mem_type,[2,4]));
nm_ucid=ucid(mem_type==0);
edges=10.^(0:0.5:4);
sust_h=[nnz(~ismember(sust_ucid,SU)),histcounts(SUC(ismember(SU,sust_ucid)),edges)]./numel(sust_ucid);
trans_h=[nnz(~ismember(trans_ucid,SU)),histcounts(SUC(ismember(SU,trans_ucid)),edges)]./numel(trans_ucid);
nm_h=[nnz(~ismember(nm_ucid,SU)),histcounts(SUC(ismember(SU,nm_ucid)),edges)]./numel(nm_ucid);
fh=figure('Color','w','Position',[100,100,240,180]);
hold on;
hs=plot(edges,sust_h,'r');
ht=plot(edges,trans_h,'b');
hn=plot(edges,nm_h,'k');
set(gca,'YScale','log','XScale','log','FontSize',10);
xlabel('In number of rings');
ylabel('Fraction of neurons');
ylim([10e-5,1]);
legend([hn,ht,hs],{'Non-memory','Transient','Sustained'});

exportgraphics(fh,fullfile('bzdata','Frac_neu_vs_num_rings_in.pdf'))

end