function nested_loop_time_const_relax_gap()
load(fullfile('binary','delay_iti_runlength_covered.mat'),'covered_tbl')
fh=figure();
hold on;
cmap=colormap('lines');
cidx=1;
covered=covered_tbl.delay;
for consecthresh=[100 200 500] % 10 tick/ms
    run_length=[];
    for jj=1:numel(covered)
        disp([consecthresh,jj])
        %pass 1
        edges = find(diff([0;covered{jj};0]==1));
        onset = edges(1:2:end-1);  % Start indices
        offset = edges(2:2:end);
        % patch through
        latmat=onset-offset.'; % row->onset, col->offset
        for onidx=1:numel(onset)
            consec=find(latmat(onidx,:)>0 & latmat(onidx,:)<=consecthresh);
            if ~isempty(consec)
                offidx=min(consec);
                covered{jj}((offset(offidx)-1):onset(onidx))=1;
            end
        end
        % pass 2
        edges = find(diff([0;covered{jj};0]==1));
        onset = edges(1:2:end-1);  % Start indices
        offset = edges(2:2:end);
        run_length =[run_length; (offset-onset)./10];  % Consecutive ones counts
    end
    chained_loops_pdf=histcounts(run_length,[0:19,20:20:180,200:100:2000],'Normalization','pdf');
    plot([0.5:19.5,30:20:190,250:100:1950],chained_loops_pdf,'-','Color',cmap(cidx,:));
    qtrs=prctile(run_length,[1,50,99]);
    xline(qtrs,'--',string(qtrs),'Color',cmap(cidx,:))
    cidx=cidx+1;
end
xlim([2,2000])
ylim([8e-6,0.1])
set(gca(),'XScale','log','YScale','log')
xlabel('Time (ms)')
ylabel('Probability density')
savefig(fh,fullfile('binary','nested_loop_time_const_relax_gap.fig'));
end
