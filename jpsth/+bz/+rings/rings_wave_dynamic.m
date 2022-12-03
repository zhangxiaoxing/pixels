load('bzdata\sums_ring_stats_all.mat')
% ring in wave

% classify rings based on su-waveid combination
% --link ring su to meta-su
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();

% wid_all=cell(1,3);
% ring_congru_idx=cell(1,3);
% for rsize=1:3
%     one_rsize=sums_all{rsize};
%     curr_sess=-1;
%     for ridx=1:size(one_rsize,1)
%         if curr_sess~=one_rsize{ridx,1}
%             curr_sess=one_rsize{ridx,1};
%             sesscid=su_meta.allcid(su_meta.sess==curr_sess);
%             sesswaveid=wrs_mux_meta.wave_id(su_meta.sess==curr_sess);
%             % update per session suid-wave map
%             sessmap=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
%         end
%         curr_waveid=cell2mat(sessmap.values(num2cell(one_rsize{ridx,3})));
%         wid_all{rsize}{ridx}=curr_waveid;
%         ring_congru_idx{rsize}{ridx}=ring_wave_type(curr_waveid);
%     end
% end

% --link ring to wave-combination congru|incongru|nonmem


% tag ring spks with trial type, trial time

% histogram of loop activity duration
mm_dur=cell(1,3);

for rsize=1:3
    one_rsize=sums_all{rsize};
    for ridx=1:size(one_rsize,1)
        rfreq_str=one_rsize{ridx,5};
        mm_dur{rsize}=[mm_dur{rsize};mean(rfreq_str.durs)./30]; % spkTS unit, i.e. 1/30 ms
    end
end

xx=5:10:95;
edges=0:10:100;
figure()
hold on
plot(xx,histcounts(mm_dur{1},edges));
plot(xx,histcounts(mm_dur{2},edges))
plot(xx,histcounts(mm_dur{3},edges))

edges=0:10:200;
per_ring_hist=struct();
[per_ring_hist.congru,per_ring_hist.incongru,per_ring_hist.nonmem,per_ring_hist.others]=deal(zeros(1,20));
for rsize=1:3
    one_rsize=sums_all{rsize};
    for ridx=1:size(one_rsize,1)
        if curr_sess~=one_rsize{ridx,1}
            curr_sess=one_rsize{ridx,1};
            sesscid=su_meta.allcid(su_meta.sess==curr_sess);
            sesswaveid=wrs_mux_meta.wave_id(su_meta.sess==curr_sess);
            sessmap=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
        end
        curr_waveid=cell2mat(sessmap.values(num2cell(one_rsize{ridx,3})));
        rwid=ring_wave_type(curr_waveid);
        rfreq_str=one_rsize{ridx,5};
        per_ring_hist.(rwid)=per_ring_hist.(rwid)+histcounts((rfreq_str.durs)./30,edges); % spkTS unit, i.e. 1/30 ms
    end
end

xx=5:10:195;
figure()
hold on
plot(xx,per_ring_hist.congru./sum(per_ring_hist.congru,'all'),'-r');
plot(xx,per_ring_hist.incongru./sum(per_ring_hist.incongru,'all'),'-b');
plot(xx,per_ring_hist.nonmem./sum(per_ring_hist.nonmem,'all'),'-k');
set(gca(),'YScale','log')
xlabel('Time (ms)')
ylabel('Probability')


function out=ring_wave_type(in)
if all(ismember(in,[1 5 7]),'all') ...
        || all(ismember(in,[2 5 8]),'all') ...
        || all(ismember(in,[3 6 7]),'all') ...
        || all(ismember(in,[4 6 8]),'all')
    out='congru';
elseif all(in==0,"all")
    out='nonmem';
elseif (any(in==5,"all") && any(in==6,"all")) ...
        ||(any(in==7,"all") && any(in==8,"all")) ...
        ||(any(ismember(in,1:2),"all") && any(in==6,"all")) ...
        ||(any(ismember(in,3:4),"all") && any(in==5,"all")) ...
        ||(any(ismember(in,[1 3]),"all") && any(in==8,"all")) ...
        ||(any(ismember(in,[2 4]),"all") && any(in==7,"all")) ...
        ||(any(in==1,"all")+any(in==2,"all")+any(in==3,"all")+any(in==4,"all")>1)
    out='incongru';
else
    out='others';
end
end
