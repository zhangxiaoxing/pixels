% [spkID,spkTS,~,~,~,~]=ephys.getSPKID_TS(18,'keep_trial',false);
% load(fullfile('binary','trials_dict.mat'),'trials_dict');
% trials=cell2mat(trials_dict(18))
function shuf_spkts=shuf_trl_ctrl(trials, spkTS, opt)
arguments
    trials
    spkTS
    opt.denovo (1,1) logical = false
end
persistent s1d3sel s1d6sel s2d3sel s2d6sel ts_cut max_t_len trials_
if opt.denovo || ~isequaln(trials_,trials)
    trials_=trials;
    processed=false(size(spkTS));
    ts_cut=nan(size(spkTS,1),2);
    max_t_len=max(diff(trials(:,1)));
    trials=[trials;trials(end,1)+max_t_len,nan(1,9)];
    for tt=1:(size(trials,1)-1)
        % disp(tt)
        curr_onset=trials(tt,1);
        next_onset=trials(tt+1,1);
        tsel=spkTS>=curr_onset & spkTS<next_onset;
        ts_cut(tsel,:)=[repmat(tt,nnz(tsel),1),spkTS(tsel)-curr_onset];
        processed(tsel)=true;
    end

    tsel=spkTS<trials(1,1);
    ts_cut(tsel,:)=[repmat(-1,nnz(tsel),1),spkTS(tsel)];

    tsel=spkTS>=trials(end,1);
    ts_cut(tsel,:)=[repmat(-2,nnz(tsel),1),spkTS(tsel)-trials(end,1)];


    s1d3sel=find(trials(:,5)==4 & trials(:,8)==3);
    s1d6sel=find(trials(:,5)==4 & trials(:,8)==6);
    s2d3sel=find(trials(:,5)==8 & trials(:,8)==3);
    s2d6sel=find(trials(:,5)==8 & trials(:,8)==6);
    ts_cut=ts_cut(isfinite(ts_cut(:,1)),:);
end
% candidate for repeats
shuf_tbl=[s1d3sel,bz.util.fisheryates(s1d3sel);...
    s1d6sel,bz.util.fisheryates(s1d6sel);...
    s2d3sel,bz.util.fisheryates(s2d3sel);...
    s2d6sel,bz.util.fisheryates(s2d6sel)];

ts_cut(:,3)=0;
for ii=reshape([s1d3sel;s1d6sel;s2d3sel;s2d6sel],1,[])
    tsel=ts_cut(:,1)==ii;
    ts_cut(tsel,3)=ts_cut(tsel,2)+trials(shuf_tbl(shuf_tbl(:,1)==ii,2),1);
end

beforets=ts_cut(ts_cut(:,1)==-1,2);
before_seg=ceil(trials(1,1)./max_t_len);
seg_len=round(trials(1,1)./before_seg);
[~,~,segid]=histcounts(beforets,0:seg_len:(before_seg*seg_len));
shuf_diff=(bz.util.fisheryates(1:before_seg).'-(1:before_seg).').*seg_len;
for ii=1:before_seg
    beforets(segid==ii)=beforets(segid==ii)+shuf_diff(ii);
end
ts_cut(ts_cut(:,1)==-1,2)=beforets;

beforets=ts_cut(ts_cut(:,1)==-1,2);
before_seg=ceil(trials(1,1)./max_t_len);
seg_len=round(trials(1,1)./before_seg);
[~,~,segid]=histcounts(beforets,0:seg_len:(before_seg*seg_len));
shuf_diff=(bz.util.fisheryates(1:before_seg).'-(1:before_seg).').*seg_len;
for ii=1:before_seg
    beforets(segid==ii)=beforets(segid==ii)+shuf_diff(ii);
end
ts_cut(ts_cut(:,1)==-1,3)=beforets;

afterts=ts_cut(ts_cut(:,1)==-2,2);
if ~isempty(afterts)
    aftermax=max(afterts);
    after_seg=ceil(aftermax./max_t_len);
    seg_len=round(aftermax./after_seg);
    [~,~,segid]=histcounts(afterts,0:seg_len:(after_seg*seg_len));
    shuf_diff=(bz.util.fisheryates(1:after_seg).'-(1:after_seg).').*seg_len;
    for ii=1:after_seg
        afterts(segid==ii)=afterts(segid==ii)+shuf_diff(ii);
    end
    ts_cut(ts_cut(:,1)==-2,3)=afterts+trials(end,1);
end
% shuf_spkts=uint32(ts_cut(:,3));
shuf_spkts=ts_cut(:,3);

end