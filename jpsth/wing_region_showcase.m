load('io_sel.mat')

wing=diff(io_entire_delay(:,[3,7]),1,2);
[~,idces]=sort(wing,'descend');
wingmat=[idces,diff(io_entire_delay(idces,[3,7]),1,2),io_entire_delay(idces,1)];
wing_reg=reg_set(idces);


load('reg_keep.mat');
greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
% any(cell2mat(cellfun(@(x) x(:,1)<100,ioselstats,'UniformOutput',false)))
nansel=io_entire_delay(:,1)<100;
sel=find(~nansel & greymatter);


rIdx=87; %87, ProS
conn_unit_all=[];
pair_unit_all=[];
in_map=[];
out_map=[];
for bin=1:6
    load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    %input
    pre_unit_sel=reg_chain_S1(:,2)==rIdx & reg_chain_S1(:,1)~=rIdx & ismember(reg_chain_S1(:,1),sel);
    pre_unit=unique(conn_chain_S1(pre_unit_sel,1));
    new_pre_unit_sel=~ismember(pre_unit,conn_unit_all);
    conn_unit_all=[conn_unit_all;pre_unit(new_pre_unit_sel)];
    
    %output
    post_unit_sel=reg_chain_S1(:,1)==rIdx & reg_chain_S1(:,2)~=rIdx & ismember(reg_chain_S1(:,2),sel);
    post_unit=unique(conn_chain_S1(post_unit_sel,2));
    new_post_unit_sel=~ismember(post_unit,conn_unit_all);
    conn_unit_all=[conn_unit_all;post_unit(new_post_unit_sel)];
    
    for uidx=1:numel(conn_unit_all)
        in_map(uidx,bin)=nnz(conn_chain_S1(pre_unit_sel,1)==conn_unit_all(uidx));
    end
    
    for uidx=1:numel(conn_unit_all)
        out_map(uidx,bin)=nnz(conn_chain_S1(post_unit_sel,2)==conn_unit_all(uidx));
    end
    
end
in_map=-abs(in_map);


connmap=[in_map,out_map];
[~,enh_idx]=sort(sum(connmap,2),'descend');
recip_sel=any(in_map(enh_idx,:),2) & any(out_map(enh_idx,:),2);
recip_idx=[enh_idx(recip_sel);enh_idx(~recip_sel)];

fh=figure('Color','w','Position',[100,100,235,235]);
ih=imagesc(connmap(recip_idx,1:6),[-8,8]);
colormap(bluewhitered(32));
colorbar()
set(gca,'XTick',[1 5])
xlabel('time (s)');
ylabel('memory neuron #');
exportgraphics(fh,'wing_input_showcase.pdf','ContentType','vector');

fh=figure('Color','w','Position',[100,100,235,235]);
ih=imagesc(connmap(recip_idx,7:12),[-8,8]);
colormap(bluewhitered(32));
colorbar()
set(gca,'XTick',[1 5])
xlabel('time (s)');
ylabel('memory neuron #');
exportgraphics(fh,'wing_output_showcase.pdf','ContentType','vector');

fh=figure('Color','w','Position',[100,100,100,235]);
hold on
per_su_sum=sum(connmap(recip_idx,:),2);
stem(find(per_su_sum>0),per_su_sum(per_su_sum>0),'Marker','none','Color','r');
stem(find(per_su_sum<0),per_su_sum(per_su_sum<0),'Marker','none','Color','b');
xlim([0.5,numel(recip_idx)+0.5]);
ylabel('net gain')
view(90,90)
set(gca,'XTick',[])
exportgraphics(fh,'wing_per_su_showcase.pdf','ContentType','vector');


keyboard;
figure('Color','w','Position',[100,100,235,235])
hold on;
% plot(-sum(in_map),'b-');
% plot(sum(out_map),'r-');
% fill([1:6,6:-1:1],[-sum(in_map),(1:6)*0],'m','FaceAlpha',0.5);
% fill([1:6,6:-1:1],[sum(out_map),fliplr(-sum(in_map))],'r','FaceAlpha',0.5);
ih=bar((1:6)-0.1,-sum(in_map),0.6,'FaceColor','b','FaceAlpha',0.5);
oh=bar((1:6)+0.1,sum(out_map),0.6,'FaceColor','r','FaceAlpha',0.5);

ylim([0,300])
xlim([0.5,6.5])
set(gca,'XTick',[1 5]);
xlabel('time bin (sec)')
ylabel('number of connection')
legend([ih,oh],{'input','output'});
exportgraphics(fh,'wing_per_bin_showcase.pdf','ContentType','vector');

return
% figure()
% hold on
% plot(input,'k-')
% plot(output,'r-')
% ylim([0,300])
