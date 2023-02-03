%% enumerate single spike and burst ring activities

% TODO: spike screen by wave / trial
function rings_burst_activiy(opt)
arguments
    opt.burst (1,1) logical = true
    opt.burstInterval (1,1) double = 600
    opt.ccg  (1,1) logical = false
    opt.append_saved (1,1) logical = false
end
% bz.rings.ring_list_bz
ppool=gcp('nocreate');
if isempty(ppool)
    if ispc
        ppool=parpool(2);
    elseif isunix
        ppool=parpool(4);
    end
end
load(fullfile('bzdata','rings_bz_wave.mat'),'rings_wave');

blame=vcs.blame();
waveids=reshape(unique(rings_wave.wave),1,[]);
sesses=reshape(unique(rings_wave.sess),1,[]);

%loop entry
if opt.append_saved
    [out,saved]=get_saved_sess();
    sesses=setdiff(sesses,saved);
else
    out=cell(0);
end
Q=parallel.FevalFuture.empty();
for sessid=sesses
    [spkID,spkTS,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
    for duration=[3 6]
        for wid=waveids
            if (contains(wid,'d3') && duration==6) ...
                    || (contains(wid,'d6') && duration ==3)
                continue
            end
            if contains(wid,'s1')
                trial_sel=find(trials(:,5)==4 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
            elseif contains(wid,'s2')
                trial_sel=find(trials(:,5)==8 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
            elseif contains(wid,'dur')
                trial_sel=find(trials(:,8)==duration & all(trials(:,9:10)>0,2));
            else
                keyboard();
            end

            sess_indices=reshape(find(rings_wave.sess==sessid & strcmp(rings_wave.wave,wid)),1,[]);
            
            for cc=sess_indices
                disp({sessid,duration,wid,cc})
                ts_id=[];
                cids=rings_wave.cids{cc};
                for in_ring_pos=1:numel(cids) 
                    one_ring_sel=spkID==cids(in_ring_pos);
                    rawts=spkTS(one_ring_sel);

                    ft_sel=strcmp(FT_SPIKE.label,num2str(cids(in_ring_pos)));
                    ft_ts=FT_SPIKE.timestamp{ft_sel};
                    ft_trl_time=FT_SPIKE.time{ft_sel};
                    ft_trl=FT_SPIKE.trial{ft_sel};

                    [~,tspos]=ismember(ft_ts,rawts);
                    ext_time=repmat(-realmax,numel(rawts),1);
                    ext_time(tspos)=ft_trl_time;

                    ext_trl=repmat(-realmax,numel(rawts),1);
                    ext_trl(tspos)=ft_trl;

                    ts_id=cat(1,ts_id,[rawts,... % 1
                        repmat(cids(in_ring_pos),numel(rawts),1),...  % 2
                        ones(numel(rawts),1)*in_ring_pos,...  % 3
                        ext_time,...  % 4
                        ext_trl]); % 5
                end
                ts_id=ts_id(ts_id(:,4)>=1 & ts_id(:,4)<(duration+1) & ismember(ts_id(:,5),trial_sel),:);
                % TODO: match trial, screen spikes
                if opt.burst
                    tsidin=cell(1,numel(cids));
                    for kk=1:numel(cids)
                        tsidin{kk}=ts_id(ts_id(:,3)==kk,1);
                    end
                    % out=relax_tag_long(in,loopIdx,recDepth,loopCnt,perSU,opt)
%                     ts=bz.rings.relax_tag_long(tsidin,[],[],[],[],"burstInterval",opt.burstInterval);
                    outkey="d"+duration+wid+"s"+sessid+"r"+numel(tsidin)+"n"+cc;
                    Q(end+1)=parfeval(ppool,@bz.rings.relax_tag_long,4,tsidin,[],[],[],[],cids,ts_id,outkey,"burstInterval",600);
%                     bz.rings.relax_tag_long(tsidin,[],[],[],[],cids,ts_id,outkey)
                else
                    error("Incomplete section")
                end
            end
        end
    end
end



%% gather from future
dbfile=fullfile("bzdata","rings_wave_burst_"+num2str(opt.burstInterval)+".db");
if ~exist(dbfile,'file')
    conn=sqlite(dbfile,"create");
    close(conn);
end
for fii=1:numel(Q)
    [ts,cids,ts_id,key]=fetchOutputs(Q(fii));
    saveOne(ts,cids,ts_id,key);
end


end
function saveOne(ts,cids,ts_id,key)
if ~isempty(ts)
    % TODO: optional remove shorter chains for each
    % onset-spike
    dbfile=fullfile("bzdata","rings_wave_burst_600.db");
    conn=sqlite(dbfile);
    for tid=1:numel(ts)
        tbl=array2table([repmat(tid,size(ts{tid},1),1),ts{tid}]);
        sqlwrite(conn,key+"_ts",tbl)
    end
    sqlwrite(conn,key+"_meta",array2table(cids))
    sqlwrite(conn,key+"_tsid",array2table(ts_id))
    close(conn);
    disp(key);
    % TODO: ccg
    %                     if opt.ccg
    %                         cursess=rings_wave.sess(ring_id);
    %                         sesspath=ephys.sessid2path(cursess);
    %                         strippath=regexp(sesspath,'(?<=\\).*','match');
    %                         sesssel=find(contains({sums_conn_str.folder},strippath));
    %                         ccg=sums_conn_str(sesssel).ccg_sc;
    %                         ccgid=sums_conn_str(sesssel).sig_con;
    %                         chainccg=[];
    %                         for ii=1:numel(cids)-1
    %                             sigsel=ccgid(:,1)==cids(ii) & ccgid(:,2)==cids(ii+1);
    %                             if nnz(sigsel)~=1
    %                                 keyboard()
    %                             end
    %                             chainccg=[chainccg;ccg(sigsel,:)];
    %                         end
    %                         out.("d"+duration).(wid).(outkey).ccgs=chainccg;
    %                     end
end
end


function [out,saved_sess]=get_saved_sess()
load('rings_wave_burst_600.mat','out');
saved_sess=[];
for dur=reshape(fieldnames(out),1,[])
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        saved_sess=[saved_sess;...
            str2double(...
            regexp(...
            fieldnames(out.(dur{1}).(wv{1})),...
            '(?<=s)\d+(?=r)','match','once'))];
    end
end
saved_sess=unique(saved_sess);
end

function plot_all(out)
figure()
hold on

perchaindur=struct();
[perchaindur.d6.size,perchaindur.d6.dur,perchaindur.d3.size,perchaindur.d3.dur,perchaindur.d6.int,perchaindur.d3.int]=deal([]);
for dur=reshape(fieldnames(out),1,[])
    %     [perchaindur.size,perchaindur.dur]=deal([]);
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            %             keyboard();
            perchaindur.(dur{1}).size=[perchaindur.(dur{1}).size,cellfun(@(x) size(x,1),out.(dur{1}).(wv{1}).(lp{1}).ts)];
            perchaindur.(dur{1}).dur=[perchaindur.(dur{1}).dur,cellfun(@(x) diff(x([1,end],3),1,1),out.(dur{1}).(wv{1}).(lp{1}).ts)./30];
            perchaindur.(dur{1}).int=[perchaindur.(dur{1}).int,cell2mat(cellfun(@(x) diff(x(:,3),1,1).',out.(dur{1}).(wv{1}).(lp{1}).ts,'UniformOutput',false))./30];
        end
    end
    statss.("d"+dur)=perchaindur;
end

d6hist=histcounts([perchaindur.d6.dur,perchaindur.d3.dur],[0:50:1000,1500,2000],'Normalization','pdf');
% d3hist=histcounts(perchaindur.d3.dur,[0:50:1000,1500,2000],'Normalization','pdf');

plot([25:50:975,1250,1750],d6hist,'k-');
% plot([25:50:975,1250,1750],d3hist,'k--');

set(gca(),'XScale','log','YScale','log');
ylim([1e-5,0.1]);
xlabel('Time (ms)')
ylabel('Probability')
end

function plot_SC(out)
end
