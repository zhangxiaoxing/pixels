% load pstats3, pstats4 from file or function
if ~exist('pstats3','var') || ~exist('pstats4','var') || ~exist('pstats5','var')
    if true
%         load('loops_proportion_stats_3.mat','pstats3');
%         load('loops_proportion_stats_4.mat','pstats4');
        load('loops_proportion_stats_all.mat','pstats3','pstats4','pstats5','mstats3','mstats4','mstats5');

    else
        if ~exist('sums_all','var')
            load(fullfile('bzdata','sums_ring_stats_all.mat'));
        end
        pstats3=load_data(3,sums_all);
        pstats4=load_data(4,sums_all);
        pstats5=load_data(5,sums_all);
        mstats3=meta_data(3,sums_all);
        mstats4=meta_data(4,sums_all);
        mstats5=meta_data(5,sums_all);
%         save('loops_proportion_stats_all.mat','mstats3','mstats4','mstats5','-append')
    end
end

sess3=cellfun(@(x) str2double(replace(x,'s','')), fieldnames(pstats3));
sess4=cellfun(@(x) str2double(replace(x,'s','')), fieldnames(pstats4));
sess5=cellfun(@(x) str2double(replace(x,'s','')), fieldnames(pstats5));
usess=unique([sess3;sess4;sess5]);
mstats=struct();
for statmtype=["congru","nonmem"]
stats=[];
if strcmp(statmtype,"congru")
    mtypeidx=1;
else
    mtypeidx=2;
end
meta=ephys.util.load_meta();
for sessid=usess.'
    sesscid=meta.allcid(meta.sess==sessid);
    sessmtype=meta.mem_type(meta.sess==sessid);
    sess_fld=sprintf('s%d',sessid);
    [spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessid);
    [ucid3,ucid4,ucid5]=deal([]);
    if isfield(pstats3,sess_fld)
        ucid3=cellfun(@(x) str2double(replace(x,'c','')),fieldnames(pstats3.(sess_fld)));
    end
    if isfield(pstats4,sess_fld)
        ucid4=cellfun(@(x) str2double(replace(x,'c','')),fieldnames(pstats4.(sess_fld)));
    end
    if isfield(pstats5,sess_fld)
        ucid5=cellfun(@(x) str2double(replace(x,'c','')),fieldnames(pstats5.(sess_fld)));
    end
    ucid=unique([ucid3;ucid4;ucid5]);
    cidmtype=arrayfun(@(x) sessmtype(sesscid==x),ucid);
    if strcmp(statmtype,"congru")
        ucid=ucid(ismember(cidmtype,1:4));
    else
        ucid=ucid(cidmtype==0);
    end
    for onecid=ucid.'
        disp([sessid,onecid]);
        su_fld=sprintf('c%d',onecid);
        totalspk=nnz(spkID==onecid);
        [max_prop,loop_spks]=deal([]);
        if isfield(pstats3,sess_fld) && isfield(pstats3.(sess_fld),su_fld)
            mtype=cell2mat(cellfun(@(x) mstats3.(sess_fld).(su_fld).(x),fieldnames(pstats3.(sess_fld).(su_fld)),'UniformOutput',false));
            if any(mtype(:,mtypeidx))
                fld=fieldnames(pstats3.(sess_fld).(su_fld));
                fldsel=fld(mtype(:,1));
                max_prop=[max_prop;cellfun(@(x) numel(pstats3.(sess_fld).(su_fld).(x)),fldsel)];
                loop_spks=[loop_spks;cell2mat(cellfun(@(x) pstats3.(sess_fld).(su_fld).(x),fldsel,'UniformOutput',false))];
            end
        end
        if isfield(pstats4,sess_fld) && isfield(pstats4.(sess_fld),su_fld)
            mtype=cell2mat(cellfun(@(x) mstats4.(sess_fld).(su_fld).(x),fieldnames(pstats4.(sess_fld).(su_fld)),'UniformOutput',false));
            if any(mtype(:,mtypeidx))
                fld=fieldnames(pstats4.(sess_fld).(su_fld));
                fldsel=fld(mtype(:,1));
                max_prop=[max_prop;cellfun(@(x) numel(pstats4.(sess_fld).(su_fld).(x)),fldsel)];
                loop_spks=[loop_spks;cell2mat(cellfun(@(x) pstats4.(sess_fld).(su_fld).(x),fldsel,'UniformOutput',false))];
            end
        end
        if isfield(pstats5,sess_fld) && isfield(pstats5.(sess_fld),su_fld)
            mtype=cell2mat(cellfun(@(x) mstats5.(sess_fld).(su_fld).(x),fieldnames(pstats5.(sess_fld).(su_fld)),'UniformOutput',false));
            if any(mtype(:,mtypeidx))
                fld=fieldnames(pstats5.(sess_fld).(su_fld));
                fldsel=fld(mtype(:,mtypeidx));
                max_prop=[max_prop;cellfun(@(x) numel(pstats5.(sess_fld).(su_fld).(x)),fldsel)];
                loop_spks=[loop_spks;cell2mat(cellfun(@(x) pstats5.(sess_fld).(su_fld).(x),fldsel,'UniformOutput',false))];
            end
        end
        if ~isempty(max_prop)
            stats=[stats;sessid,onecid,totalspk,numel(max_prop),max(max_prop),numel(unique(loop_spks))];
        end
    end
end

stats(:,7)=stats(:,5)./stats(:,3).*100;
stats(:,8)=stats(:,6)./stats(:,3).*100;
if strcmp(statmtype,"congru")
    mstats.congru=stats;
else
    mstats.nonmem=stats;
end
end
save('loops_proportion_stats.mat','mstats')

% singleci=bootci(500,@(x) mean(x),stats(:,7));
% alterci=bootci(500,@(x) mean(x),stats(:,8));
% mm=mean(stats(:,7:8));

perbin=[];
anovamat=[];
for bin=20:20:60
    cbinsel=mstats.congru(:,4)>bin-20 & mstats.congru(:,4)<=bin;
    cmm=mean(mstats.congru(cbinsel,8));
    cci=bootci(500,@(x) mean(x),mstats.congru(cbinsel,8));
    cc=mstats.congru(cbinsel,8);
    cc(:,2)=bin;
    cc(:,3)=1;
    
    nbinsel=mstats.nonmem(:,4)>bin-20 & mstats.nonmem(:,4)<=bin;
    nmm=mean(mstats.nonmem(nbinsel,8));
    nci=bootci(500,@(x) mean(x),mstats.nonmem(nbinsel,8));
    nn=mstats.nonmem(nbinsel,8);
    nn(:,2)=bin;
    nn(:,3)=0;
    
    anovamat=[anovamat;cc;nn];
    perbin=[perbin;mean(cbinsel),cmm,cci(1),cci(2),mean(nbinsel),nmm,nci(1),nci(2)];
end
fh=figure('Color','w','Position',[32,32,160,225]);
hold on
ch=plot(10:20:60,perbin(:,2),'-r');
fill([10:20:60,fliplr(10:20:60)],[perbin(:,3);flip(perbin(:,4))],'r','EdgeColor','none','FaceAlpha',0.1)

nh=plot(10:20:60,perbin(:,6),'-k');
fill([10:20:60,fliplr(10:20:60)],[perbin(:,7);flip(perbin(:,8))],'k','EdgeColor','none','FaceAlpha',0.1)

ylabel('Loops-associated spikes (%)')
xlabel('Number of loops involved')
legend([ch,nh],{'Same-memory','Non-memory'},'Location','northoutside');
set(gca(),'XTick',0:20:60,'YTick',0:10:max(ylim()))
xlim([0,60])
exportgraphics(fh,'loops_spk_proportion.pdf');
anovan(anovamat(:,1),{anovamat(:,2),anovamat(:,3)})

% 
% fh=figure('Color','w','Position',[32,32,195,235]);
% hold on;
% swarmchart([ones(size(stats,1),1);2*ones(size(stats,1),1)],[stats(:,7);stats(:,8)],1,'o','filled','MarkerFaceColor','r','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none','XJitterWidth',0.5)
% % bar(1:2,mm,0.6,'FaceColor','none','EdgeColor','k','LineWidth',1);
% errorbar(1:2,mm,[singleci(1),alterci(1)]-mm,[singleci(2),alterci(2)]-mm,'k.','CapSize',15);
% ylabel('Loops-associated spikes (%)');
% set(gca(),'XTick',1:2,'XTickLabel',{'Single loop','Alternative paths'},'XTickLabelRotation',90)
% xlim([0.5,2.5]);
% ylim([0,60]);
% exportgraphics(fh,'loops_spk_proportion.pdf');

function pstats=load_data(rsize,sums_all)
% if rsize==4
%     load(fullfile('bzdata','sums_ring_stats_4.mat'),'sums4');
%     rstats=sums4;
% else
%     load(fullfile('bzdata','sums_ring_stats_3.mat'),'sums3');
%     rstats=sums3;
% end
rstats=sums_all{rsize-2};
rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==rsize,:);

usess=unique(cell2mat(rstats(:,1)));
pstats=struct();
for sessid=usess.'
    pstats.(sprintf('s%d',sessid))=struct();
    [spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessid);
    rid=find(cell2mat(rstats(:,1))==sessid);
    for ri=reshape(rid,1,[])
        disp([sessid,ri]);
        ts_id=[];
        cids=rstats{ri,3};
        per_cid_spk_cnt=cids;
        for in_ring_pos=1:numel(cids) % TODO, 1:rsize
            one_ring_sel=spkID==cids(in_ring_pos);
            per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
            ts_id=cat(1,ts_id,[spkTS(one_ring_sel),ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
        end
        ts_id=sortrows(ts_id,1);
        ts_id=[ts_id,rstats{ri,5}.tags];

        for cidx=1:numel(rstats{ri,3})
            cid=rstats{ri,3}(cidx);
            if ~isfield(pstats.(sprintf('s%d',sessid)),sprintf('c%d',cid))
                pstats.(sprintf('s%d',sessid)).(sprintf('c%d',cid))=struct();
            end
            pstats.(sprintf('s%d',sessid)).(sprintf('c%d',cid)).(sprintf('r%d',ri))=ts_id(ts_id(:,2)==cidx & ts_id(:,3)~=0,1);
        end
    end
end
end

function mstats=meta_data(rsize,sums_all)
meta=ephys.util.load_meta();
rstats=sums_all{rsize-2};
rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==rsize,:);

usess=unique(cell2mat(rstats(:,1)));
mstats=struct();
for sessid=usess.'
    sesscid=meta.allcid(meta.sess==sessid);
    sessmtype=meta.mem_type(meta.sess==sessid);
    mstats.(sprintf('s%d',sessid))=struct();
    rid=find(cell2mat(rstats(:,1))==sessid);
    for ri=reshape(rid,1,[])
        disp([sessid,ri]);
        %find mem type here
        rtypes=arrayfun(@(x) sessmtype(sesscid==x),rstats{ri,3});
        congru=all(ismember(rtypes,1:2),'all') | all(ismember(rtypes,3:4),'all');
        nonmem=all(rtypes==0,'all');
        for cidx=1:numel(rstats{ri,3})
            cid=rstats{ri,3}(cidx);
            if ~isfield(mstats.(sprintf('s%d',sessid)),sprintf('c%d',cid))
                mstats.(sprintf('s%d',sessid)).(sprintf('c%d',cid))=struct();
            end
            mstats.(sprintf('s%d',sessid)).(sprintf('c%d',cid)).(sprintf('r%d',ri))=[congru,nonmem];
        end
    end
end
end

