keyboard()
bz=false;
if false
    %both data result from su_region_recip.m
    load('rings.mat')
%     load('candidate_count.mat')
end
if true
    load('loops_count.mat')
else
    shufrpt=size(rings_shuf,1);
%     sesscnt=countConnedSession();
    for midx=1:size(rings,1)
        for sess=1:size(rings,2)
            for bin=1:size(rings,3)
                motif_count(midx,sess,bin)=size(cell2mat([(rings(midx,sess,bin,1)),(rings(midx,sess,bin,2))]'),1);
            end
        end
    end

    for midx=1:size(rings,1)
        for sess=1:size(rings,2)
            for bin=1:size(rings,3)
                motif_inact_count(midx,sess,bin)=size(cell2mat([(rings_inact(midx,sess,bin,1)),(rings_inact(midx,sess,bin,2))]'),1);
            end
        end
    end

    motif_count_shuf=nan(size(rings,1),size(rings,2),size(rings,3),shufrpt);
    for midx=1:size(rings,1)
        for sess=1:size(rings,2)
            for bin=1:size(rings,3)
                for rpt=1:shufrpt
                    motif_count_shuf(midx,sess,bin,rpt)=size(cell2mat([squeeze(rings_shuf(rpt,midx,sess,bin,1));squeeze(rings_shuf(rpt,midx,sess,bin,2))]),1);
                end
            end
        end
    end
    
    
    motif_inact_count_shuf=nan(size(rings,1),size(rings,2),size(rings,3),shufrpt);
    for midx=1:size(rings,1)
        for sess=1:size(rings,2)
            for bin=1:size(rings,3)
                for rpt=1:shufrpt
                    motif_inact_count_shuf(midx,sess,bin,rpt)=size(cell2mat([squeeze(rings_shuf_inact(rpt,midx,sess,bin,1));squeeze(rings_shuf_inact(rpt,midx,sess,bin,2))]),1);
                end
            end
        end
    end
    
    %% baseline
    if ~bz
        for midx=1:size(base_rings,1)
            for sess=1:size(base_rings,2)
                base_motif_count(midx,sess)=size(cell2mat([(base_rings(midx,sess,1)),(base_rings(midx,sess,2))]'),1);
            end
        end

        base_motif_count_shuf=nan(size(base_rings,1),size(base_rings,2),shufrpt);
        for midx=1:size(base_rings,1)
            for sess=1:size(base_rings,2)
                for rpt=1:shufrpt
                    base_motif_count_shuf(midx,sess,rpt)=size(cell2mat([squeeze(base_rings_shuf(rpt,midx,sess,1));squeeze(base_rings_shuf(rpt,midx,sess,2))]),1);
                end
            end
        end
        % shuf_mm=mean(shuf_count,3);
        %% get figure
    else
        base_motif_count=[];
        base_motif_count_shuf;
    end
    save('loops_count.mat','motif_count','motif_count_shuf','motif_inact_count','motif_inact_count_shuf','base_motif_count','base_motif_count_shuf')
end
    close all
for msize=3:5
%     plotOne(msize,motif_inact_count,motif_count_shuf,base_motif_count,base_motif_count_shuf,true);
    plotOne(msize,motif_count,motif_count_shuf,base_motif_count,base_motif_count_shuf,false);
end

return

function plotOne(msize,motif_delay, motif_shuf,motif_base,motif_base_shuf,inact)
midx=msize-2;
mmshufC=squeeze(sum(motif_shuf(midx,:,:,:),2));
shufci=bootci(1000,@(x) mean(x),mmshufC');
shufmm=mean(mmshufC,2);
% mmbase_shuf=mmbase_shuf_C./sum(candidate_count_base(:,midx));
% mmshuf=mmshufC./sum(candidate_delay(:,:,midx));


fh=figure('Color','w','Position',[100,100,150,150]);
hold on;
if inact 
    
    shufcib=bootci(1000,@(x) mean(x),squeeze(sum(motif_base_shuf(midx,:,:),2)));
    ci(:,1)=sum(motif_base(midx,:),2)./shufcib;
    mm=[sum(motif_base(midx,:),2)./mean(size(sum(motif_base_shuf(midx,:,:),2)));nan;squeeze(sum(motif_delay(midx,:,:),2))./shufmm];
    ci(:,3:8)=squeeze(sum(motif_delay(midx,:,:),2))'./shufci;
    rh=bar(-1:6,mm,0.6,'FaceColor','k','FaceAlpha',0.2);
    errorbar(-1:6,mm,ci(1,:)-mm',ci(2,:)-mm','k.','CapSize',3,'MarkerSize',2);
    ylim([0.5,5000])
else
    
    ci(:,3:8)=squeeze(sum(motif_delay(midx,:,:),2))'./shufci;
    mm=[nan;nan;squeeze(sum(motif_delay(midx,:,:),2))./shufmm];
    rh=bar(-1:6,mm,0.6,'FaceColor','w','EdgeColor','k');
    errorbar(-1:6,mm,ci(1,:)-mm',ci(2,:)-mm','k.','CapSize',3,'MarkerSize',2);
    ylim([0.5,5000])
end

disp([nanmean(mm),nanstd(mm)/sqrt(nnz(~isnan(mm)))])

xlim([-1.75,6.75])
yline(1,'--','Color',[0.5,0.5,0.5])
set(gca,'XTick',[-1,0 5],'XTickLabel',{'ITI','0','5'},'YScale','log','YColor','k','FontSize',10,'Ytick',[1,10,100,1000]);
xlabel('Time (sec)')
ylabel('Observed / expected')
% keyboard()
if inact
    exportgraphics(fh,sprintf('motif_of_%d_inact.pdf',msize),'ContentType','vector');
else
    exportgraphics(fh,sprintf('motif_of_%d.pdf',msize),'ContentType','vector');
end

if ~evalin('base','bz')
    hold on
    if inact
        return
        [pref_rings,nonp_rings]=collect_data(midx,true);
        fp=errorbar([-1,1:6],mean(pref_rings(:,[1,4:9])./msize),std(pref_rings(:,[1,4:9])./msize)./sqrt(size(pref_rings,1)),'-ro','MarkerSize',3);
        fn=errorbar([-1,1:6],mean(nonp_rings(:,[1,4:9])./msize),std(nonp_rings(:,[1,4:9])./msize)./sqrt(size(nonp_rings,1)),'-bo','MarkerSize',3);
%         fp=plot([-1,1:6],mean(pref_rings(:,[1,4:9])),'-ro','MarkerSize',3);
%         fn=plot([-1,1:6],mean(nonp_rings(:,[1,4:9])),':ro','MarkerSize',3);
    else
        fh=figure('Color','w','Position',[100,100,150,150]);
        hold on;
        [pref_rings,nonp_rings]=collect_data(midx);
        fp=errorbar(1:6,nanmean(pref_rings(:,4:9)./msize),nanstd(pref_rings(:,4:9)./msize)./sqrt(sum(~isnan(pref_rings(:,4:9)))),'-ro','MarkerSize',3);
        fn=errorbar(1:6,nanmean(nonp_rings(:,4:9)./msize),nanstd(nonp_rings(:,4:9)./msize)./sqrt(sum(~isnan(nonp_rings(:,4:9)))),'-bo','MarkerSize',3);
%         plot(1:6,mean(pref_rings(:,4:9)),'-ro','MarkerSize',3);
%         plot(1:6,mean(nonp_rings(:,4:9)),':ro','MarkerSize',3);
    end
    
    if inact
        ylim([0,1]);
    else
        ylim([0,2]);
    end
    set(gca,'YColor','r','FontSize',10,'YTick',0:max(ylim()));
    ylabel('Ring cycle frequency')
%     legend([fp,fn],{'Prefered','Non-pref.'})
    text(0,0.5,sprintf('n=%d',size(pref_rings,1)),'Color','r');
    set(gca,'XTick',[0,5],'XTickLabel',[0,5],'YScale','linear','YColor','k','FontSize',10,'Ytick',0:0.5:max(ylim()));
%     ylim([0,0.5])
    xlim([0,6.5])
    xlabel('Time (sec)')
end
% legend([rh,sh],{'recorded','shuffled'});
% keyboard()
if inact
    exportgraphics(fh,sprintf('motif_of_%d_inact_freq.pdf',msize),'ContentType','vector');
else
    exportgraphics(fh,sprintf('motif_of_%d_freq.pdf',msize),'ContentType','vector');
end
end


function sesscnt=countConnedSession()
fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0114_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
sesscnt=false(114,1);
for I=1:114
    disp(I)
    lbound=100000*I;
    ubound=100000*(I+1);
    for bin=1:6
        sel11=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
        %                 keyboard()
        if nnz(sel11)>0
            sesscnt(I)=true;
        end
    end
end
end


function [pref_rings,nonp_rings,neu_count]=collect_data(midx,inact)
datapath='ring_freq';
if exist('inact','var') && inact
    flist=ls(fullfile(datapath,sprintf('*freq_inact_%d_*.mat',midx+2)));
else
    flist=ls(fullfile(datapath,sprintf('*freq_%d_*.mat',midx+2)));
end

ring_list=cell(0,8);
for i=1:size(flist,1)
    fn=fullfile(datapath,deblank(flist(i,:)));
    
    if ~isfile(fn)
        continue
    end
    load(fn,'all_rings');
    ring_list=[ring_list;all_rings];
end

% load significant3_ringlist.mat
% r=sort(reshape(cell2mat(ring_list(:,3)),3,[])'+cell2mat(ring_list(:,2))*100000,2);
% shuf_sel=ismember(r,significant3_ringlist,'rows');
% ring_list=ring_list(shuf_sel,:);

% ringOccur=sum(cell2mat(ring_list(:,5:6)),2);
% totalTrial=sum(cellfun(@(x) size(x,1), ring_list(:,7:8)),2);
% ring_list=ring_list(ringOccur>=totalTrial,:);

occur_sel=(([ring_list{:,4}]==1 & [ring_list{:,5}]>cellfun(@(x) size(x,1),ring_list(:,7))') |...
([ring_list{:,4}]==2 & [ring_list{:,6}]>cellfun(@(x) size(x,1),ring_list(:,8))'));
ring_list=ring_list(occur_sel',:);

neu_count=numel(unique(reshape(cell2mat(ring_list(:,3)),numel(ring_list{1,3}),[])'+cell2mat(ring_list(:,2))*100000));
pref_rings=[];
nonp_rings=[];
for i=1:size(ring_list,1)
    if ring_list{i,4}==1
        prefhist=mean(ring_list{i,7});
        nonphist=mean(ring_list{i,8});
    else
        prefhist=mean(ring_list{i,8});
        nonphist=mean(ring_list{i,7});
    end
    pref_rings=[pref_rings;prefhist];
    nonp_rings=[nonp_rings;nonphist];    
end

end



