%idces from hebb_pattern_showcase.m
if ~exist('rings','var')
    addpath(fullfile('npy-matlab-master','npy-matlab'))
    addpath('fieldtrip-20200320')
    ft_defaults
    load rings.mat
    load 114_sorted_file_path.mat
    delay=6;
    if isunix
        cd('~/pixels/jpsth')
        homedir='/home/zx/neupix/wyt';
    elseif ispc
        homedir='k:\neupix\wyt';
    end
end


%midx=3;
r1=[];
for bb=1:6
    r1=[r1;cell2mat(rings(midx,:,bb,1)')];
end
r2=[];
for bb=1:6
    r2=[r2;cell2mat(rings(midx,:,bb,2)')];
end
r1(:,midx+3)=1;r2(:,midx+3)=2;ringsm=[r1(:,1:(midx+2));r2(:,1:(midx+2))];samps=[r1(:,midx+3);r2(:,midx+3)];
keyboard

for ssidx=srange
    if ssidx>length(ringsm)
        break
    end
    sessIdx=idivide(int32(ringsm(ssidx,1)),int32(100000));
    sample=samps(ssidx);


    folder=regexp(sorted_fpath{sessIdx},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
    currmodel='selec';
    sustIds=[];nonselIds=[];
    transIds=int32(rem(ringsm(ssidx,:),100000))';
    
    %% %%%%%% localize %%%%%%%%%%%%%
    freg=load('0831_selec_conn_chain_duo_6s_1_2.mat','pair_reg','pair_chain');
    load reg_keep.mat
    suid=ringsm(ssidx,:);
    reg=[];
    for i1=suid(:)'
        reg=[reg;{i1,reg_set{freg.pair_reg(find(freg.pair_chain==i1,1))}}];
    end
    %% end of localize
    
    msize=numel(transIds);
    [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
    if sample==1
        cfg.trials = find(spktrial.trialinfo(:,5)==4 & spktrial.trialinfo(:,8)==delay);
    else
        cfg.trials = find(spktrial.trialinfo(:,5)==8 & spktrial.trialinfo(:,8)==delay);
    end
    examples=cell(0);
    cocount=[];
    for tidx=cfg.trials(:)'
        disp(tidx);
        ts_id=[];
        for seqid=1:msize
            tsel=spktrial.trial{seqid}==tidx;
            ts=spktrial.time{seqid}(tsel);
            ts=ts(ts>1 & ts<7);
            ts(2,:)=seqid;
            ts_id=[ts_id;ts'];
        end
        [~,s]=sort(ts_id(:,1));
        ts_id=ts_id(s,:);
        ts_id_tagged=relax_tag(ts_id,msize);
        coact_count=sum(ts_id_tagged(:,3));
        if coact_count>100
            examples{end+1}=ts_id_tagged;
            cocount(end+1,:)=[coact_count,tidx];
        end
    end
    if ~isempty(cocount)
        [~,s]=sort(cocount(:,1),'descend');
        for i=1
            plotseq(examples{s(i)},msize,cocount(s(i),2),transIds,sessIdx,ssidx,cocount(s(i),1));
        end
    end
end


function out=tag_ring(in,msize)
out=in;
out(:,3)=0;
if msize==3
    for i=1:length(in)-3
        if (isequal(in(i:i+3,2),[1 2 3 1]') || isequal(in(i:i+3,2),[2 3 1 2]') || isequal(in(i:i+3,2),[3 1 2 3]')) ...
                && max(diff(in(i:i+3,1)))<0.01 && min(diff(in(i:i+3,1)))>0.001
            out(i:i+3,3)=1;
        end
    end
elseif msize==4
    for i=1:length(in)-4
        if (isequal(in(i:i+4,2),[1 2 3 4 1]')...
                || isequal(in(i:i+4,2),[2 3 4 1 2]')...
                || isequal(in(i:i+4,2),[3 4 1 2 3]')...
                || isequal(in(i:i+4,2),[4 1 2 3 4]'))...
                && max(diff(in(i:i+4,1)))<0.01 ...
                && min(diff(in(i:i+4,1)))>0.001
            out(i:i+3,3)=1;
        end
    end
end
end



function plotseq(in,msize,tidx,ids,sessIdx,ssidx,cocount)
fh=figure('Color','w','Position',[100,100,400,400]);
hold on
color={'r','b','c','m','g'};
ph=matlab.graphics.chart.primitive.Line.empty(0,5);
for bin=1:6
    for i=1:msize
        j=i;
        selp=in(:,2)==i & in(:,3)==1 & in(:,1)>=bin & in(:,1)<bin+1;
        seln=in(:,2)==i & in(:,3)==0 & in(:,1)>=bin & in(:,1)<bin+1;
        tph=plot(repmat(in(selp,1)',2,1)-bin,(bin-1)*(msize+1)+repmat([j-0.6;j+0.6],1,nnz(selp)),'-','Color',color{i});
%         ph(i)=tph(1);
        plot(repmat(in(seln,1)',2,1)-bin,(bin-1)*(msize+1)+repmat([j-0.6;j+0.6],1,nnz(seln)),'-','Color',[0.8,0.8,0.8]);

    end
end
xlabel('Time within time-bin(s)')
ylabel('Time-bin')
set(gca,'XTick',0:0.5:1,'YTick',(msize/2):(msize+1):((msize+1)*6),'YTickLabel',1:6)
if msize==3
    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d',sessIdx,tidx,ids(1),ids(2),ids(3)));
    %     print(fh,sprintf('spike_seq_%d_%d_%d.png',ids(1),ids(2),ids(3)),'-dpng');
elseif msize==4
    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d, %d',sessIdx,tidx,ids(1),ids(2),ids(3),ids(4)));
    %     print(fh,sprintf('spike_seq_%d_%d_%d_%d.png',ids(1),ids(2),ids(3),ids(4)),'-dpng');
elseif msize==4
    title(sprintf('Sess#%d, Trial#%d',sessIdx,tidx));    
end
print(fh,sprintf('r%d_spike_seq_%03d_%d.png',msize,cocount,ssidx),'-dpng');
if false
%     legend(ph,arrayfun(@(x) ['cell ',num2str(x)],1:5,'UniformOutput',false),'NumColumns',3,'Location','northoutside')
    ylim([0,6*(msize+1)]);
    fh.Position(3:4)=[225,225];
    set(gca,'FontSize',10);
    exportgraphics(fh,sprintf('r%d_spike_seq_%d_%d.pdf',msize,cocount,ssidx),'ContentType','vector');
end

close(fh);
end

