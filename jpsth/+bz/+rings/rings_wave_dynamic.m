% ring in wave
% load('bzdata\sums_ring_stats_all.mat')
function  rings_wave_dynamic(sums_all)
% TODO optional global filter to remove spikes outside delay

% classify rings based on su-waveid combination
% --link ring su to meta-su
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();

% --link ring to wave-combination congru|incongru|nonmem

% tag ring spks with trial type, trial time

% histogram of loop activity duration
if false
    mm_dur=[];

    for rsize=1:3
        one_rsize=sums_all{rsize};
        for ridx=1:size(one_rsize,1)
            rfreq_str=one_rsize{ridx,5};
            mm_dur=[mm_dur;mean(rfreq_str.durs)./30]; % spkTS unit, i.e. 1/30 ms
        end
    end

    xx=5:10:195;
    edges=0:10:200;
    figure()
    hold on
    plot(xx,histcounts(mm_dur,edges));
    set(gca(),'YScale','log')
    xlabel('Time (ms)')
    ylabel('Probability')


    return

    mm_dur=cell(1,3);

    for rsize=1:3
        one_rsize=sums_all{rsize};
        for ridx=1:size(one_rsize,1)
            rfreq_str=one_rsize{ridx,5};
            mm_dur{rsize}=[mm_dur{rsize};mean(rfreq_str.durs)./30]; % spkTS unit, i.e. 1/30 ms
        end
    end

    xx=5:10:195;
    edges=0:10:200;
    figure()
    hold on
    h3=plot(xx,histcounts(mm_dur{1},edges),'k');
    h4=plot(xx,histcounts(mm_dur{2},edges),'b');
    h5=plot(xx,histcounts(mm_dur{3},edges),'r');
    legend([h3,h4,h5],{'3-Neuron','4-Neuron','5-Neuron'})
    set(gca(),'YScale','log')
end
curr_sess=-1;
edges=[0:1:10,20:10:200];
per_ring_hist=struct();
[per_ring_hist.congru,per_ring_hist.incongru,per_ring_hist.nonmem,per_ring_hist.others]=deal(zeros(1,29));
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
        rwid=bz.rings.ring_wave_type(curr_waveid);
        rfreq_str=one_rsize{ridx,5};
        per_ring_hist.(rwid)=per_ring_hist.(rwid)+histcounts((rfreq_str.durs)./30,edges); % spkTS unit, i.e. 1/30 ms
    end
end

%%
load('sums_conn.mat','sums_conn_str');
qcmat=cell2mat({sums_conn_str.qc}.'); % from bz.sums_conn -> bz.goodccg
fwhm_all=(qcmat(:,2)-250)./30; % offset left half of symmatric ccg
fcxx=0:0.8:20;
disp('F.C. mean, sd, sem')
disp([mean(fwhm_all), std(fwhm_all), std(fwhm_all)./sqrt(numel(fwhm_all))])
fcyy=histcounts(fwhm_all,fcxx);
fciqr=prctile(fwhm_all,[25,50,75]);
%%

xx=[0.5:1:9.5,15:10:195];
figure()
hold on
% plot(xx,per_ring_hist.congru./sum(per_ring_hist.congru,'all'),'-r');
% plot(xx,per_ring_hist.incongru./sum(per_ring_hist.incongru,'all'),'-b');
% plot(xx,per_ring_hist.nonmem./sum(per_ring_hist.nonmem,'all'),'-k');
count_sum=per_ring_hist.congru+per_ring_hist.incongru+per_ring_hist.nonmem+per_ring_hist.others;
plot(xx,count_sum./sum(count_sum,'all'),'-r');
plot(0.4:0.8:19.6,fcyy./sum(fcyy,'all'),'-k');
%
% xline(fciqr(2),'k-')
% xline(fciqr(2),'k--')
set(gca(),'YScale','log','XScale','log')
xlim([0.3,200])
ylim([1e-3,1])
xlabel('Time (ms)')
ylabel('Probability')
end
