% tag chain-like spike sequences in population recording data
% sust(ain)=burst
% TODO: rename

% load(fullfile('bzdata','chains_mix.mat'),'chains_uf');
% chains=chains_uf;
% clear chains_uf;

% TODO: wave consistant, wave inconsistent

function out=chain_sust_tag(chains,opt)
arguments
    chains
    opt.ccg (1,1) logical = true
    opt.burstInterval (1,1) double = 300
end

% tsidin={0,25,50,75,100,125,150};
% ts=chain_recur(tsidin,[],[],1);
% keyboard()

if opt.ccg
    load('sums_conn_10.mat','sums_conn_str');
end
%% set up
% all_chains=fieldnames(pstats.congru);
ch_len=cellfun(@(x) numel(x),chains.cids);
waveids=reshape(unique(chains.wave(ch_len>4)),1,[]);
sesses=reshape(unique(chains.sess(ch_len>4)),1,[]);

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

            sess_indices=reshape(find(chains.sess==sessid & strcmp(chains.wave,wid) & ch_len>4),1,[]);

            for cc=sess_indices
                ts_id=[];
                cids=chains.cids{cc};
                %                 disp({sessid,ri});
                for in_chain_pos=1:numel(cids) 
                    one_chain_sel=spkID==cids(in_chain_pos);
                    rawts=spkTS(one_chain_sel);

                    ft_sel=strcmp(FT_SPIKE.label,num2str(cids(in_chain_pos)));
                    ft_ts=FT_SPIKE.timestamp{ft_sel};
                    ft_trl_time=FT_SPIKE.time{ft_sel};
                    ft_trl=FT_SPIKE.trial{ft_sel};

                    [~,tspos]=ismember(ft_ts,rawts);
                    ext_time=repmat(-realmax,numel(rawts),1);
                    ext_time(tspos)=ft_trl_time;

                    ext_trl=repmat(-realmax,numel(rawts),1);
                    ext_trl(tspos)=ft_trl;

                    ts_id=cat(1,ts_id,[rawts,... % 1
                        repmat(cids(in_chain_pos),numel(rawts),1),...  % 2
                        ones(numel(rawts),1)*in_chain_pos,...  % 3
                        ext_time,...  % 4
                        ext_trl]); % 5
                end
                % optional remove non-wave spikes
                ts_id=ts_id(ts_id(:,4)>=1 & ts_id(:,4)<(duration+1) & ismember(ts_id(:,5),trial_sel),:);
                %                 ts_id=sortrows(ts_id,1);
                tsidin=cell(0);
                for kk=1:numel(cids)
                    tsidin{kk}=ts_id(ts_id(:,3)==kk,1);
                end
                ts=chain_recur(tsidin,[],[],1,'burstInterval',opt.burstInterval);
                if ~isempty(ts)

                    % TODO: optional remove shorter chains for each
                    % onset-spike

                    outkey="s"+sessid+"c"+cc;
                    out.("d"+duration).(wid).(outkey).ts=ts;
                    out.("d"+duration).(wid).(outkey).meta={cids,chains.tcoms(cc)};
                    out.("d"+duration).(wid).(outkey).ts_id=ts_id;

                    % TODO: ccg
                    if opt.ccg
                        cursess=chains.sess(cc);
                        sesspath=ephys.sessid2path(cursess);
                        strippath=regexp(sesspath,'(?<=\\).*','match');
                        sesssel=find(contains({sums_conn_str.folder},strippath));
                        ccg=sums_conn_str(sesssel).ccg_sc;
                        ccgid=sums_conn_str(sesssel).sig_con;
                        chainccg=[];
                        for ii=1:numel(cids)-1
                            sigsel=ccgid(:,1)==cids(ii) & ccgid(:,2)==cids(ii+1);
                            if nnz(sigsel)~=1
                                keyboard()
                            end
                            chainccg=[chainccg;ccg(sigsel,:)];
                        end
                        out.("d"+duration).(wid).(outkey).ccgs=chainccg;
                    end
                end
            end
        end
    end
end


for dur=reshape(fieldnames(out),1,[])
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            ccgs=out.(dur{1}).(wv{1}).(lp{1}).ccgs;
            ccg_qual=[];
            for ii=1:size(ccgs,1)
                ccg_qual=[ccg_qual;bz.good_ccg(ccgs(ii,:))];
                %1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising
                %edge 6:falling edge
            end
            out.(dur{1}).(wv{1}).(lp{1}).ccg_qual=ccg_qual;
        end
    end
end

blame=vcs.blame();
save(fullfile("bzdata","chain_sust_tag_"+opt.burstInterval+".mat"),'out','blame','-v7.3')
stats(out);
end

function stats(out)
% load(fullfile("bzdata","chain_sust_tag_"+burstInterval+".mat"),'out')
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

d6hist=histcounts([perchaindur.d3.dur,perchaindur.d6.dur],[0:20:100,200:100:600],'Normalization','pdf');
% d3hist=histcounts(perchaindur.d3.dur,0:50:600,'Normalization','pdf');
figure()
hold on
plot([10:20:90,150:100:550],d6hist,'-k');
% plot(25:50:575,d3hist,'--k');
set(gca(),'XScale','log','YScale','log');




end


function out=chain_recur(in,chainDepth,currIdx,recDepth,opt)
arguments
    in
    chainDepth
    currIdx
    recDepth
    opt.burstInterval (1,1) double = 300 % ticks, at 30k sps
end
% disp("recDepth"+num2str(recDepth))
% disp("interval "+num2str(opt.burstInterval));
out=cell(0); %{[depth,id]}
clen=size(in,2);

if isempty(chainDepth)
    chainDepth=1;
    currIdx=1;
end

while currIdx<=numel(in{chainDepth})
    % same neuron
    if currIdx<=numel(in{chainDepth})-1 && in{chainDepth}(currIdx+1)-in{chainDepth}(currIdx)<opt.burstInterval
        link=[chainDepth,currIdx,in{chainDepth}(currIdx)];
        rec=chain_recur(in,chainDepth,currIdx+1,recDepth+1,'burstInterval',opt.burstInterval);
%         if chainDepth==1 && currIdx==1
%             keyboard()
%         end
        for ii=1:numel(rec)
            out(end+1)={[link;rec{ii}]};
        end
    else % no branch last ts
        if chainDepth==clen
            out={[chainDepth,currIdx,in{chainDepth}(currIdx)]};
        end
    end
    % branch
    if chainDepth<clen
        nxtd=chainDepth+1;
        nidx=1;
        while nidx<=numel(in{nxtd}) && in{nxtd}(nidx)<in{chainDepth}(currIdx)+24
            nidx=nidx+1;
        end
        if nidx<=numel(in{nxtd}) && in{nxtd}(nidx)<in{chainDepth}(currIdx)+300 % in window
            link=[chainDepth,currIdx,in{chainDepth}(currIdx)];
            rec=chain_recur(in,nxtd,nidx,recDepth+1,'burstInterval',opt.burstInterval);
            for ii=1:numel(rec)
                out(end+1)={[link;rec{ii}]};
            end
        end
    end
    if recDepth==1
        while currIdx<=numel(in{1})-1 && in{1}(currIdx+1)<in{1}(currIdx)+opt.burstInterval
            currIdx=currIdx+1;
        end
        currIdx=currIdx+1;
    else
        break
    end
end
end

function plot_all()
figure()
hold on
intvs=[150 300 600];
colors=['k','b'];
for pidx=1:3
    load(fullfile("bzdata","chain_sust_tag_"+intvs(pidx)+".mat"),'out')
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

    d6hist=histcounts(perchaindur.d6.dur,[0:50:1000,1500,2000],'Normalization','pdf');
    d3hist=histcounts(perchaindur.d3.dur,[0:50:1000,1500,2000],'Normalization','pdf');
    plot([25:50:975,1250,1750],d6hist,'-','Color',colors(pidx));
    plot([25:50:975,1250,1750],d3hist,'--','Color',colors(pidx));
end
set(gca(),'XScale','log','YScale','log');
ylim([1e-4,0.1]);
xlim([0,1000]);
xlabel('Time (ms)')
ylabel('Probability')
end

