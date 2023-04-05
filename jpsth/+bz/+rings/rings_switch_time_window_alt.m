% Primary developing branch since 2023.1.2
function [single_su_multi_ring,ssmr_meta]=rings_switch_time_window_alt(single_su_multi_ring)
%data file from rings_time_constant.m
if false
    load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats')
    rstats.congru=pstats.congru;
    pstats=rstats;
    clear rstats
end
if false
    % for each session
    single_su_multi_ring.d3=struct();
    single_su_multi_ring.d6=struct();
    ssmr_meta.d3=struct();
    ssmr_meta.d6=struct();

    all_rings=fieldnames(pstats.congru);
    waveids={[1 5],[2 5],[3 6],[4 6],[1 7],[2 8],[3 7],[4 8],[7 8]};
    sess=str2double(string(regexp(all_rings,'(?<=s)\d+(?=r.*)','match')));
    for wid=1:8
        waveid=waveids{wid};
        for sessid=reshape(unique(sess),1,[])
            sess_ring=all_rings(sess==sessid);
            trials=pstats.congru.(sess_ring{1}).trials;
            for duration=[3 6]
                if any(ismember([1 2 5],waveid),'all')
                    trial_sel=find(trials(:,5)==4 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
                elseif any(ismember([3 4 6],waveid),'all')
                    trial_sel=find(trials(:,5)==8 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
                elseif all(ismember(waveid,7:8),"all")
                    trial_sel=find(trials(:,8)==duration & all(trials(:,9:10)>0,2));
                else
                    keyboard();
                end

                for rr=reshape(sess_ring,1,[])
                    if ~all(ismember(pstats.congru.(rr{1}).rstats{4},waveid),2) ...
                            || (duration==3 && any(ismember(pstats.congru.(rr{1}).rstats{4},[2,4,8]),"all")) ...
                            || (duration==6 && any(ismember(pstats.congru.(rr{1}).rstats{4},[1,3,7]),"all"))
                        continue
                    end
                    ts_id=pstats.congru.(rr{1}).ts_id;
                    rcids=pstats.congru.(rr{1}).rstats{3};
                    for su=rcids
                        sutag="s"+sessid+"w"+wid+"u"+num2str(su);
                        if ~isfield(single_su_multi_ring.("d"+num2str(duration)),sutag)
                            single_su_multi_ring.("d"+num2str(duration)).(sutag)=ts_id(ts_id(:,2)==su & ismember(ts_id(:,5),trial_sel),4:5);
                            ssmr_meta.("d"+num2str(duration)).(sutag)={pstats.congru.(rr{1}).rstats};
                        end
                        single_su_multi_ring.("d"+num2str(duration)).(sutag)(:,end+1)=ts_id(ts_id(:,2)==su & ismember(ts_id(:,5),trial_sel) ,6);
                        ssmr_meta.("d"+num2str(duration)).(sutag)(end+1)={pstats.congru.(rr{1}).rstats};
                    end
                end
            end
        end
    end

    blame=vcs.blame();
    save(fullfile("bzdata","single_su_multi_ring.mat"),"ssmr_meta","single_su_multi_ring","blame");
else
    load(fullfile("bzdata","single_su_multi_ring.mat"),"ssmr_meta","single_su_multi_ring");
end
if false
    out=stats(single_su_multi_ring);
    assignin('base','out',out)
    blame=vcs.blame();
    save('rings_switch_time_window.mat','out','blame')
else
    load('rings_switch_time_window.mat','out')
end
sel1=out.act_cnt_1(:,1)>0;
data1=out.act_cnt_1(sel1,:)./out.bin_cnt_1(sel1,:);
mm1=mean(data1).*100;
sem1=std(data1)./sqrt(size(data1,1)).*100;

sel10=out.act_cnt_10(:,1)>0;
data10=out.act_cnt_10(sel10,:)./out.bin_cnt_10(sel10,:);
mm10=mean(data10).*100;
sem10=std(data10)./sqrt(size(data1,10)).*100;

sel50=out.act_cnt_50(:,1)>0;
data50=out.act_cnt_50(sel50,:)./out.bin_cnt_50(sel50,:);
mm50=mean(data50).*100;
sem50=std(data50)./sqrt(size(data1,50)).*100;

figure()
hold on
pxx=0.025:0.05:0.975;
fill([pxx,fliplr(pxx)],[mm1-sem1,fliplr(mm1+sem1)],'k','EdgeColor','none','FaceAlpha',0.2);
fill([pxx,fliplr(pxx)],[mm10-sem10,fliplr(mm10+sem10)],'b','EdgeColor','none','FaceAlpha',0.2);
fill([pxx,fliplr(pxx)],[mm50-sem50,fliplr(mm50+sem50)],'r','EdgeColor','none','FaceAlpha',0.2);
h1=plot(pxx,mm1,'k-');
h10=plot(pxx,mm10,'b-');
h50=plot(pxx,mm50,'r-');
legend([h1,h10,h50],{'Single communal loop','10 communal loops','50 communal loops'},'Location','northoutside','Orientation','horizontal');
xlabel('Time (s)')
ylabel('Probability of activation by communal loops (%)');
set(gca(),'XTick',0:0.25:1,'YTick',0:25:100)

end

function out=stats(single_su_multi_ring)
% TODO showcase

% TODO histogram
% 3s and 6s
cnt3s=cellfun(@(x) size(x,2), struct2cell(single_su_multi_ring.d3));
cnt6s=cellfun(@(x) size(x,2), struct2cell(single_su_multi_ring.d6));
% hist3=histcounts(cnt3s,[0:10:100,500],'Normalization','pdf');
hhist=histcounts([cnt3s;cnt6s],[0.5:200.5],'Normalization','cdf');
figure()
hold on;
% plot(5:10:95,mean([hist3(1:end-1);hist6(1:end-1)]),'-k');
% plot(95:10:105,mean([hist3(end-1:end);hist6(end-1:end)]),'--k');
plot([1:200],hhist,'-k');
% set(gca(),'YScale','log','XTick',0:50:150)
% ylim([4e-4,0.1])
xlim([0,175])
xlabel('Number of composite loops')
ylabel('Probability density')
title('Single SU in composite loops')


% TODO prob of switch vs interval
[out.bin_cnt_1,out.act_cnt_1,out.bin_cnt_10,out.act_cnt_10,out.bin_cnt_50,out.act_cnt_50]=deal([]);
for dur=[3 6]
    durtag="d"+num2str(dur);
    for onesu=reshape(fieldnames(single_su_multi_ring.(durtag)),1,[])
        disp(onesu{1})
        sudata=single_su_multi_ring.(durtag).(onesu{1});

%         total_len=numel(unique(sudata(:,2)))*dur; % total time in second
%         total_len./sum(sudata(:,3:end)>0);

        sudata(:,1)=sudata(:,1)-1;
%         keyboard()
        [bin_cnt_1,act_cnt_1,bin_cnt_10,act_cnt_10,bin_cnt_50,act_cnt_50]=prob_switch_loop(sudata,dur);
        out.bin_cnt_1=[out.bin_cnt_1;bin_cnt_1];
        out.act_cnt_1=[out.act_cnt_1;act_cnt_1];
        out.bin_cnt_10=[out.bin_cnt_10;bin_cnt_10];
        out.act_cnt_10=[out.act_cnt_10;act_cnt_10];
        out.bin_cnt_50=[out.bin_cnt_50;bin_cnt_50];
        out.act_cnt_50=[out.act_cnt_50;act_cnt_50];
    end
end

end

function [bin_cnt_1,act_cnt_1,bin_cnt_10,act_cnt_10,bin_cnt_50,act_cnt_50]=prob_switch_loop(sudata,dur)
delta=1/30000;
win_wid=1;
step=0.05;
bin_cnt_1=zeros(1,win_wid/step);
act_cnt_1=zeros(1,win_wid/step);

bin_cnt_10=zeros(1,win_wid/step);
act_cnt_10=zeros(1,win_wid/step);

bin_cnt_50=zeros(1,win_wid/step);
act_cnt_50=zeros(1,win_wid/step);

trls=reshape(unique(sudata(:,2)),1,[]);
for trl=trls
    for lead=3:size(sudata,2)
        if size(sudata,2)>13 && rem(lead,10)==0
            disp([trl,lead,size(sudata,2)])
        end

        lead_act=sudata(sudata(:,2)==trl,lead);
        if nnz(lead_act)==0
            continue
        end
        all_act=unique(lead_act(lead_act>0));
        for one_act=reshape(all_act,1,[])
            one_time=max(sudata(sudata(:,lead)==one_act,1));
            if one_time+delta>dur-1
                continue
            end
            for follow=setdiff(3:size(sudata,2),lead)
                follow_time=sudata(sudata(:,2)==trl & sudata(:,follow)~=0,1);
                bin_edge=(one_time+delta):step:(one_time+win_wid+delta);
                one_cnt=histcounts(follow_time,bin_edge,"Normalization","cumcount")>0;
%                 if nnz(one_cnt)>2
%                     keyboard()
%                 end
                bin_cnt_1=bin_cnt_1+1;
                act_cnt_1=act_cnt_1+one_cnt;
            end
            % 10 loop
            if size(sudata,2)>13
                follow_set=setdiff(3:size(sudata,2),lead);
                for rpt=1:10
                    one_rpt_set=randsample(follow_set,10);
                    follow_time=sudata(sudata(:,2)==trl & any(sudata(:,one_rpt_set)~=0,2),1);
                    bin_edge=(one_time+delta):step:(one_time+win_wid+delta);
                    one_cnt=histcounts(follow_time,bin_edge,"Normalization","cumcount")>0;
                    bin_cnt_10=bin_cnt_10+1;
                    act_cnt_10=act_cnt_10+one_cnt;
                end
            end
            % 50 loop
            if size(sudata,2)>53
                follow_set=setdiff(3:size(sudata,2),lead);
                for rpt=1:10
                    one_rpt_set=randsample(follow_set,50);
                    follow_time=sudata(sudata(:,2)==trl & any(sudata(:,one_rpt_set)~=0,2),1);
                    bin_edge=(one_time+delta):step:(one_time+win_wid+delta);
                    one_cnt=histcounts(follow_time,bin_edge,"Normalization","cumcount")>0;
                    bin_cnt_50=bin_cnt_50+1;
                    act_cnt_50=act_cnt_50+one_cnt;
                end
            end

        end
    end
end

end


