% motifs PSTH

sessid=100;
maxid=4;
[~,ttrial]=wave.chain_loop_SC_spk(false,sessid,maxid,'skip_plot',true);


covered=false(1,1);
tsmin=max(min(cell2mat(ttrial(:,4))));
figure()
tiledlayout(2,1)
% motif raster
nexttile()
hold on
% for each motif
for mmid=1:size(ttrial,1)
%extract ts
ts=(ttrial{mmid,4}(:,2)-tsmin)./30;
%scatter %optional alternative plot
scatter(ts,mmid,'k','|')
onset=floor(ts(1))+1;
offset=ceil(ts(end));
covered(onset:offset)=true;
end

xlim([0,ceil(numel(covered)./1000)*1000]);
xlabel('Time (s)')
ylim([0.5,size(ttrial,1)+5.5])
ylabel('Motifs #')
edges = find(diff([0,covered,0]==1));
onset = edges(1:2:end-1);  % Start indices
run_length = edges(2:2:end)-onset;  % Consecutive ones counts
ymax=max(ylim());
for ii=1:numel(onset)
    plot([onset(ii),onset(ii)+run_length(ii)],[ymax,ymax],'r-','LineWidth',4);
end
% PSTH
nexttile()
hold on