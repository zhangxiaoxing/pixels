% ring in wave
% load('bzdata\sums_ring_stats_all.mat')
function  rings_wave_dynamic(sums_all)
% TODO optional global filter to remove spikes outside delay
global_init();
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

%% individual loops
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

%% individual FC
load('sums_conn.mat','sums_conn_str');
qcmat=cell2mat({sums_conn_str.qc}.'); % from bz.sums_conn -> bz.goodccg
fwhm_all=(qcmat(:,2)-250)./30; % offset left half of symmatric ccg
fcxx=0:0.8:20;
disp('F.C. mean, sd, sem')
disp([mean(fwhm_all), std(fwhm_all), std(fwhm_all)./sqrt(numel(fwhm_all))])
fcyy=histcounts(fwhm_all,fcxx);
fciqr=prctile(fwhm_all,[25,50,75]);

%% wave bump width
if false
    bump3=[];
    bump6=[];
    com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false);
    for fn=reshape(fieldnames(com_map),1,[])
        for subfn=reshape(fieldnames(com_map.(fn{1})),1,[])
            if isfield(com_map.(fn{1}).(subfn{1}),'fwhm3')
                bump3=[bump3,cell2mat(com_map.(fn{1}).(subfn{1}).fwhm3.values())];
            end
            if isfield(com_map.(fn{1}).(subfn{1}),'fwhm6')
                bump6=[bump6,cell2mat(com_map.(fn{1}).(subfn{1}).fwhm6.values())];
            end
        end
    end

    bxx=0.125:0.25:6.25;
    bump3hist=histcounts(bump3,bxx,'Normalization','probability');
    bump6hist=histcounts(bump6,bxx,'Normalization','probability');
end
%% chains end-to-end span
load('chains_mix.mat','chains_uf');
chains=chains_uf;
clear chains_uf;
wids=reshape(unique(chains.wave),1,[]);
t_span=struct();
for dur=3:3:6
    for wid=wids
        lsel=cellfun(@(x) numel(x)>4,chains.cids);
        wsel=chains.wave==wid & chains.dur==dur;
        dtcom=cellfun(@(x) x(end)-x(1),chains.tcoms(lsel & wsel));
%         if ~isfield(t_span,wid+num2str(dur))
%             t_span.(wid+num2str(dur))=[];
%         end
        if ~isempty(dtcom)
            t_span.("d"+num2str(dur)).(wid)=dtcom;
        end
    end
end

cxx=[0:200:1000,1500:500:6000];
chist3=histcounts([t_span.d3.olf_s1;t_span.d3.olf_s2;t_span.d3.s1d3;t_span.d3.s2d3]*1000/4,cxx,'Normalization','probability');
chist6=histcounts([t_span.d6.olf_s1;t_span.d6.olf_s2;t_span.d6.s1d6;t_span.d6.s2d6]*1000/4,cxx,'Normalization','probability');


%% composite loops

hrstats=load(fullfile('bzdata','hebbian_ring.mat'),'stats');
C=struct2cell(hrstats.stats).';
expd=[C{:}];
expdd=[expd{:}];
hrhist=histcounts(expdd,0:10:250,'Normalization','probability');


%%
xx=[0.5:1:9.5,15:10:195];
figure()
hold on
% plot(xx,per_ring_hist.congru./sum(per_ring_hist.congru,'all'),'-r');
% plot(xx,per_ring_hist.incongru./sum(per_ring_hist.incongru,'all'),'-b');
% plot(xx,per_ring_hist.nonmem./sum(per_ring_hist.nonmem,'all'),'-k');
count_sum=per_ring_hist.congru+per_ring_hist.incongru+per_ring_hist.nonmem+per_ring_hist.others;
if false
    fch=plot(0.4:0.8:19.6,fcyy./sum(fcyy,'all'),'--k.');
    looph=plot(xx,count_sum./sum(count_sum,'all'),'-kx');
    hrh=plot(5:10:245,hrhist,'-.k+');
    if false
        b3h=plot((0.25:0.25:6).*1000,bump3hist,'-b');
        b6h=plot((0.25:0.25:6).*1000,bump6hist,'-c');
    end
    p_cxx=[100:200:900,1250:500:5750];%cxx=[0:200:2000,2500:500:6000];
    c3h=plot(p_cxx,chist3,'-k^');
    c6h=plot(p_cxx,chist6,'--ko');
else
    fch=plot(0.4:0.8:19.6,fcyy./sum(fcyy,'all'),'-','Color',"#D35400");
    looph=plot(xx,count_sum./sum(count_sum,'all'),'-','Color',"#27AE60");
    hrh=plot(5:10:245,hrhist,'-','Color',"#8E44AD");
    if false
        b3h=plot((0.25:0.25:6).*1000,bump3hist,'-b');
        b6h=plot((0.25:0.25:6).*1000,bump6hist,'-c');
    end
    p_cxx=[100:200:900,1250:500:5750];%cxx=[0:200:2000,2500:500:6000];
    c3h=plot(p_cxx,chist3,'-','Color',"#2980B9");
    c6h=plot(p_cxx,chist6,'-','Color',"#C0392B");
end
xline(3000,'--k')
xline(6000,'--k')
%
% xline(fciqr(2),'k-')
% xline(fciqr(2),'k--')
set(gca(),'YScale','log','XScale','log')
xlim([0.3,6000])
ylim([1e-3,1])
xlabel('Time (ms)')
ylabel('Probability')
legend([fch,looph,hrh,c3h,c6h],{'FC','Loop','Composite loops','3s chain-loop','6s chain-loop'},'Location','northoutside','Orientation','horizontal')
end
