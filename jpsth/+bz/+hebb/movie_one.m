%sess,activ3,tags3,all3_sel,activ4,tags4,all4_sel,selidx,ii
% bz.hebb.movie_one(sess,activ3,tags3,all3_sel,activ4,tags4,all4_sel,selidx,ii);
function movie_one(sess,activ3,tags3,all3_sel,activ4,tags4,all4_sel,selidx,hbb)
fhh=figure('Color','w');
% for idx=0:6
%     updateHebb(fhh,idx)
%     pause(0.5)
% end

[spkID,spkTS,trials,~,~]=ephys.getSPKID_TS(sess);
ring_spk_sums=cell(0);
suseq=[];
%% ring of 3
for rr=reshape(find(all3_sel),1,[])
    ts_id=[];
    for iid=activ3(rr,:) % TODO, 1:rsize
        spk=spkTS(spkID==iid);
        spk(:,2)=iid;
        ts_id=[ts_id;spk];
    end
    ring_spk_sums=[ring_spk_sums;sortrows(ts_id,1),full(tags3{rr}.tags)];
    suseq=[suseq;activ3(rr,:),nan];
end
%% ring of 4
for rr=reshape(find(all4_sel),1,[])
    ts_id=[];
    for iid=activ4(rr,:) % TODO, 1:rsize
        spk=spkTS(spkID==iid);
        spk(:,2)=iid;
        ts_id=[ts_id;spk];
    end
    ring_spk_sums=[ring_spk_sums;sortrows(ts_id,1),full(tags4{rr}.tags)];
    suseq=[suseq;activ4(rr,:)];
end

trial_sel=find(trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0 & trials(:,5)==selidx*4);

for tt=141%reshape(trial_sel,1,[])
    tonset=trials(tt,1)+30000;
    toffset=trials(tt,2);
    
    for onset=2.57*30000+tonset%((0:0.05:5.9)*30000+tonset)

%         onset=tonset+5.05*30000;
        offset=onset+0.075*30000;
        per_ring_cnt=cellfun(@(x) nnz(x(x(:,1)>=onset & x(:,1)<offset,3)>0),ring_spk_sums);
%         per_ring_t_cnt=cellfun(@(x) nnz(x(x(:,1)>=tonset & x(:,1)<toffset,3)>0),ring_spk_sums);
        if nnz(per_ring_cnt)<3 || sum(per_ring_cnt)<11
            continue
        end

        colors=distinguishable_colors(nnz(per_ring_cnt));
        %%
        fsh=figure('Color','w','Position',[10,100,350,200]);
        hold on
        cidx=1;
        ringtag=[];
        % for neuron view
        suids=reshape(unique(suseq(isfinite(suseq))),1,[]);
        spkinrings=containers.Map('KeyType','double','ValueType','any');
        spkalls=containers.Map('KeyType','double','ValueType','any');

        for jj=suids
            spkinrings(jj)=[];
        end

        for ii=1:numel(per_ring_cnt)
            if per_ring_cnt(ii)==0, continue;end
            rcolor=colors(cidx,:);
            cidx=cidx+1;
            ts_id=ring_spk_sums{ii};
            tsel=ts_id(:,1)>=onset & ts_id(:,1)<offset;
            % for rings view
            seq=suseq(ii,isfinite(suseq(ii,:)));
            for jj=seq
                xin=ts_id(tsel & ts_id(:,2)==jj & ts_id(:,3)>0,1);
                spkinrings(jj)=[spkinrings(jj);xin];
            end

        end

        yidx=1;
        for jj=[375 189 475 535 524 171]+10000%suids
            spkall=spkTS(spkID==jj & spkTS>=onset & spkTS<offset);
            
            xin=spkinrings(jj);
            spkalls(yidx)=[jj;xin];
            xout=spkall(~ismember(spkall,xin));

            if rem(yidx,2)==0
                fill([onset,offset,offset,onset],[-0.5,-0.5,0.5,0.5]+yidx,'k','EdgeColor','none','FaceAlpha',0.05);
            end

            plot(repmat(xin,1,2).',repmat([yidx-0.5;yidx+0.5],1,numel(xin)),'-k','LineWidth',1);
            yidx=yidx+1;
        end

%         xlim([tonset+5.05*30000,tonset+5.2*30000])
%         set(gca(),'XTick',tonset:1500:toffset+15000,'XTickLabel',0:50:6000,'YTick',[]);
        set(gca(),'XTick',onset:1500:offset,'XTickLabel',0:50:100,'YTick',[]);
        xlabel('Time (ms)')

        title(sprintf('Sess%d, HebbPatt%d, Trial%d, spk%d',sess,hbb,tt,sum(per_ring_cnt)))

        hists=cell2mat(cellfun(@(x) histcounts(x(2:end),onset:6:offset),spkalls.values(),'UniformOutput',false).');
        sp=0;
        figure(fhh)
        v=VideoWriter('hebb_sc.mp4','MPEG-4');
        open(v);
        updateHebb(fhh,sp);
        box off
        set(gca,'Visible',false)
        th=text(250,640,'0 msec','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',20);
        writeVideo(v,getframe(fhh))
        for ii=1:size(hists,2)
            if any(hists(:,ii))
                sp=sp+1;
                updateHebb(fhh,sp);
            end
            box off
        set(gca,'Visible',false)
        delete(th)
        th=text(250,640,sprintf('%0.1f msec',ii*0.2),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',20);
        writeVideo(v,getframe(fhh))
        end
        close(v)
    end
end


    


end
function updateHebb(fh,idx)
persistent imdata
if isempty(imdata)
    imdata=cell(15,1);
    for ii=0:14
        imdata{ii+1}=imread(fullfile('SC','hebb',sprintf('hebb%d.png',ii)));
    end
end
figure(fh)
image(imdata{idx+1});
truesize(fh)
end

