function fr_sel_idx(opt)
arguments
    opt.delay (1,1) double {mustBeMember(opt.delay,[3 6])} =6
    opt.per_bin (1,1) logical = false
end
meta=ephys.util.load_meta();
homedir=ephys.util.getHomedir('type','raw');

[selidx,eselidx]=deal([]);

for sess=1:116
    sesspath=meta.allpath{find(meta.sess==sess,1)};
    fpath=fullfile(homedir,sesspath,'FR_All_1000.hdf5');
    fpath=replace(fpath,'\',filesep());
    fr=h5read(fpath,'/FR_All');
    trial=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    s1sel=(trial(:,5)==4 & trial(:,8)==opt.delay & trial(:,9)>0 & trial(:,10)>0);
    s2sel=(trial(:,5)==8 & trial(:,8)==opt.delay & trial(:,9)>0 & trial(:,10)>0);
    
    %% each neuron
    s1uids=meta.allcid(meta.sess==sess & ismember(meta.mem_type,1:2).');
    [~,loc]=ismember(s1uids,suid);
    if opt.per_bin
        bins=meta.per_bin(:,meta.sess==sess & ismember(meta.mem_type,1:2).');
        s1us1=arrayfun(@(x) mean(fr(s1sel,loc(x),find(bins(:,x))+4),'all') ,1:numel(loc));
        s1us2=arrayfun(@(x) mean(fr(s2sel,loc(x),find(bins(:,x))+4),'all') ,1:numel(loc));
    else
        s1us1=mean(fr(s1sel,loc,5:10),[1 3]);
        s1us2=mean(fr(s2sel,loc,5:10),[1 3]);
    end
    selidx=[selidx,(s1us1-s1us2)./(s1us1+s1us2)];

    s2uids=meta.allcid(meta.sess==sess & ismember(meta.mem_type,3:4).');
    [~,loc]=ismember(s2uids,suid);
    if opt.per_bin
        bins=meta.per_bin(:,meta.sess==sess & ismember(meta.mem_type,3:4).');
        s2us1=arrayfun(@(x) mean(fr(s1sel,loc(x),find(bins(:,x))+4),'all') ,1:numel(loc));
        s2us2=arrayfun(@(x) mean(fr(s2sel,loc(x),find(bins(:,x))+4),'all') ,1:numel(loc));        
    else
        s2us1=mean(fr(s1sel,loc,5:10),[1 3]);
        s2us2=mean(fr(s2sel,loc,5:10),[1 3]);
    end
    selidx=[selidx,(s2us2-s2us1)./(s2us1+s2us2)];

    e1sel=(trial(:,5)==4 & trial(:,8)==opt.delay & trial(:,10)==0);
    e2sel=(trial(:,5)==8 & trial(:,8)==opt.delay & trial(:,10)==0);
    
    if nnz(e1sel)<5 || nnz(e2sel)<5
        continue;
    end

    if opt.per_bin
        [~,loc]=ismember(s1uids,suid);
        bins=meta.per_bin(:,meta.sess==sess & ismember(meta.mem_type,1:2).');
        s1ue1=arrayfun(@(x) mean(fr(e1sel,loc(x),find(bins(:,x))+4),'all') ,1:numel(loc));
        s1ue2=arrayfun(@(x) mean(fr(e2sel,loc(x),find(bins(:,x))+4),'all') ,1:numel(loc));
    else
        s1ue1=mean(fr(e1sel,loc,5:10),[1 3]);
        s1ue2=mean(fr(e2sel,loc,5:10),[1 3]);
    end
    interm_idx=(s1ue1-s1ue2)./(s1ue1+s1ue2);
    interm_idx((s1ue1+s1ue2)==0)=0;
    eselidx=[eselidx,interm_idx];

    if opt.per_bin
        [~,loc]=ismember(s2uids,suid);
        bins=meta.per_bin(:,meta.sess==sess & ismember(meta.mem_type,3:4).');
        s2ue1=arrayfun(@(x) mean(fr(e1sel,loc(x),find(bins(:,x))+4),'all') ,1:numel(loc));
        s2ue2=arrayfun(@(x) mean(fr(e2sel,loc(x),find(bins(:,x))+4),'all') ,1:numel(loc));
    else
        s2ue1=mean(fr(e1sel,loc,5:10),[1 3]);
        s2ue2=mean(fr(e2sel,loc,5:10),[1 3]);
    end
    interm_idx=(s2ue2-s2ue1)./(s2ue2+s2ue1);
    interm_idx((s2ue1+s2ue2)==0)=0;
    eselidx=[eselidx,interm_idx];

end
sem=[std(selidx)./sqrt(numel(selidx)),std(eselidx)./sqrt(numel(eselidx))];
mm=[mean(selidx),mean(eselidx)];
fh=figure('Color','w','Position',[32,32,230,230]);
% swarmchart(ones(size(selidx)),selidx,1,'k','filled','o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
hold on
bar(1:2,mm,'FaceColor','w','EdgeColor','k');
errorbar(1:2,mm,sem,'k.','CapSize',15);
set(gca(),'XTick',1:2, 'XTickLabel',{'Correct','Error'})
ylabel('Selectivity index')
if opt.per_bin
    ylim([0,0.3])
    exportgraphics(fh,'fr_sel_idx_per_bin.pdf')
else
    ylim([-0.02,0.12])
    exportgraphics(fh,'fr_sel_idx.pdf')
end

end


