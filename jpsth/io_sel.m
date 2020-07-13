load('reg_keep.mat');
bin=1;
load(sprintf('0712_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
in_out_sel=nan(length(reg_set),12);
for reg_idx=1:length(reg_set)
    pair_fw= nnz((pair_reg(:,1)~=reg_idx) & (pair_reg(:,2)==reg_idx));
    pair_rev= nnz((pair_reg(:,1)==reg_idx) & (pair_reg(:,2)~=reg_idx));
    pair_count=pair_fw+pair_rev;
    
    in_conn_S1= nnz((reg_chain_S1(:,1)~=reg_idx) & (reg_chain_S1(:,2)==reg_idx));
    in_sel_S1=nnz((reg_chain_S1(:,1)~=reg_idx) & (reg_chain_S1(:,2)==reg_idx) & pref_chain_S1(:,bin)>0);

    out_conn_S1= nnz((reg_chain_S1(:,2)~=reg_idx) & (reg_chain_S1(:,1)==reg_idx));
    out_sel_S1=nnz((reg_chain_S1(:,2)~=reg_idx) & (reg_chain_S1(:,1)==reg_idx) & pref_chain_S1(:,bin+6)>0);
    
    auto_pair=nnz((pair_reg(:,1)==reg_idx) & (pair_reg(:,2)==reg_idx))*2;
    auto_conn_S1= nnz((reg_chain_S1(:,1)==reg_idx) & (reg_chain_S1(:,2)==reg_idx));
    
    in_out_sel(reg_idx,:)=[pair_count,in_conn_S1,in_conn_S1/pair_count, ...%1 2 3
        in_sel_S1,in_sel_S1/pair_count,...% 4 5
        out_conn_S1,out_conn_S1/pair_count,...% 6 7
        out_sel_S1,out_sel_S1/pair_count,...% 8 9
        auto_pair,auto_conn_S1,auto_conn_S1/auto_pair]; % 10 11 12
end

countsel=in_out_sel(:,1)>100;
iosel=in_out_sel(countsel,:);
sink_sel=false(size(iosel,1),1);
source_sel=false(size(iosel,1),1);
for i=1:size(iosel,1)
    currcount=iosel(i,1);
    [~,~,p]=crosstab(1:2*currcount>currcount,[1:currcount>iosel(i,4),1:currcount>iosel(i,8)]);
%     chisqsel(i)=p*numel(chisqsel)<0.05; %bonferroni adjust
    if p<0.05 && iosel(i,8)>iosel(i,4)
        source_sel(i)=true;
    elseif p<0.05 && iosel(i,8)<iosel(i,4)
        sink_sel(i)=true;
    end
end


regsel=reg_set(countsel);
figure('Color','w');
hold on
scatter(iosel(~(sink_sel|source_sel),5),iosel(~(sink_sel|source_sel),9),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
scatter(iosel(source_sel,5),iosel(source_sel,9),'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.5);
scatter(iosel(sink_sel,5),iosel(sink_sel,9),'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha',0.5);
disp('source')
disp(regsel(source_sel));
disp('sink')
disp(regsel(sink_sel));
plot([0,1],[0,1],'k:');
xlim([0,0.15]);
ylim([0,0.15]);
xlabel('in-density');
ylabel('out-density');