load('io_sel.mat'); %save('io_sel.mat','ioselstats','io_entire_delay','io_early_delay','io_late_delay','reg_set');
%%%%%%%%%%%%%%%%%NEW CODE 0724%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% abundened due to distractions in STR-CTX function conn
% regline=textread('K:\neupix\meta\regClass.csv','%s');
% regclass=cell(0);
% regclassId=[];
% for i=1:length(regline)
%     regclassT=regexp(regline{i},'(.*),(?>STR|TH|CTX|ELSE)','tokens','once');
%     if endsWith(regline{i},',CTX')
%         regclass{end+1}=regclassT;
%         regclassId(end+1)=0;
%     elseif endsWith(regline{i},',TH')
%         regclass{end+1}=regclassT;
%         regclassId(end+1)=1;
%     elseif endsWith(regline{i},',STR')
%         regclass{end+1}=regclassT;
%         regclassId(end+1)=2;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        

%%%%%reminder of data tag
%     in_out_sel(reg_idx,:)=[pair_count,in_conn_S1,in_conn_S1/pair_count, ...%1 2 3
%         in_sel_S1,in_sel_S1/pair_count,...% 4 5
%         out_conn_S1,out_conn_S1/pair_count,...% 6 7
%         out_sel_S1,out_sel_S1/pair_count,...% 8 9
%         auto_pair,auto_conn_S1,auto_conn_S1/auto_pair]; % 10 11 12

ctdPath='K:\code\ctd\correct_error\showcase';
ctdf=fullfile(ctdPath,'corr_error_decoding.hdf5');
c_reg_raw=deblank(h5read(ctdf,'/c_reg'));
e_reg=deblank(h5read(ctdf,'/e_reg'));
c_mat_raw=h5read(ctdf,'/c_mat');
e_mat=h5read(ctdf,'/e_mat');

sel=ismember(c_reg_raw,e_reg);
c_reg=c_reg_raw(sel);
c_mat=c_mat_raw(:,:,:,sel);

% all(strcmp(c_reg,e_reg)) %confirmed to be true
dec_accu=nan(length(c_reg),6);
corr_prefer_dec=[];
corr_conn_dec=[];
corr_auto_prefer=[];
corr_auto_conn=[];
corr_auto_dec=[];
corr_reg=cell(0);


for reg=1:length(c_reg)
    mmc=nonzeros(mean(c_mat(:,:,:,reg),3).*eye(32))./100;
    mme=nonzeros(mean(e_mat(:,:,:,reg),3).*eye(32))./100;
    dec_accu(reg,:)=[mean(mmc(9:32)),mean(mmc(9:20)),mean(mmc(21:32)),...
        mean(mme(9:32)),mean(mme(9:20)),mean(mme(21:32))];
    % delay corr, early corr, late corr, delay err, early err, late err
    
    io_reg_idx=find(strcmp(reg_set,c_reg{reg}));
    if ~isempty(io_reg_idx)
        corr_reg{end+1}=c_reg{reg};
        corr_prefer_dec(end+1,:)=[diff(dec_accu(reg,[4 1]),1,2),diff(io_entire_delay(io_reg_idx,[5,9]),1,2),...
            diff(dec_accu(reg,[5 2]),1,2),diff(io_early_delay(io_reg_idx,[5,9]),1,2),...
            diff(dec_accu(reg,[6 3]),1,2),diff(io_late_delay(io_reg_idx,[5,9]),1,2)];
        corr_conn_dec(end+1,:)=[diff(dec_accu(reg,[4 1]),1,2),diff(io_entire_delay(io_reg_idx,[3,7]),1,2),...
            diff(dec_accu(reg,[5 2]),1,2),diff(io_early_delay(io_reg_idx,[3,7]),1,2),...
            diff(dec_accu(reg,[6 3]),1,2),diff(io_late_delay(io_reg_idx,[3,7]),1,2)];
        corr_auto_prefer(end+1,:)=[];
        
        % ctd-currDecodingIdx,io_diff_idx
    end
end

fh=figure('Color','w','Position',[100,100,280,280]);
% subplot(1,3,1)
sh=scatter(corr_prefer_dec(:,1),corr_prefer_dec(:,2),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
for i=1:length(corr_reg)
    text(corr_prefer_dec(i,1),corr_prefer_dec(i,2),corr_reg{i},'HorizontalAlignment','center')
end
[r,p]=corrcoef(corr_prefer_dec(:,1),corr_prefer_dec(:,2));
legend(sh,sprintf('r=%.3f, p=%.3f',r(1,2),p(1,2)));
xlabel('delta decoding accuracy')
ylabel('delta input-output selectivity')
title('entire delay');
keyboard
exportgraphics(fh,'io_sel_decoding_corr.pdf','ContentType','vector');

fh=figure('Color','w','Position',[100,100,280,280]);
% subplot(1,3,1)
sh=scatter(corr_conn_dec(:,1),corr_conn_dec(:,2),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
for i=1:length(corr_reg)
    text(corr_conn_dec(i,1),corr_conn_dec(i,2),corr_reg{i},'HorizontalAlignment','center')
end
[r,p]=corrcoef(corr_conn_dec(:,1),corr_conn_dec(:,2));
legend(sh,sprintf('r=%.3f, p=%.3f',r(1,2),p(1,2)));
xlabel('delta decoding accuracy')
ylabel('delta input-output selectivity')
title('entire delay');
keyboard
 exportgraphics(fh,'io_conn_decoding_corr.pdf','ContentType','vector');




% print('io_selectivity_decoding_corr.png','-dpng','-r300');

figure('Color','w','Position',[100,100,800,400]);
subplot(1,2,1)
sh=scatter(corr_prefer_dec(:,3),corr_prefer_dec(:,4),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
[r,p]=corrcoef(corr_prefer_dec(:,3),corr_prefer_dec(:,4));
legend([sh],sprintf('r=%.3f, p=%.3f',r(1,2),p(1,2)));
xlabel('delta decoding accuracy')
ylabel('delta input-output selectivity')
title('early delay');

subplot(1,2,2)
sh=scatter(corr_prefer_dec(:,5),corr_prefer_dec(:,6),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
[r,p]=corrcoef(corr_prefer_dec(:,5),corr_prefer_dec(:,6));
legend([sh],sprintf('r=%.3f, p=%.3f',r(1,2),p(1,2)));
xlabel('delta decoding accuracy')
ylabel('delta input-output selectivity')
title('late delay');

print('io_selectivity_decoding_corr.png','-dpng','-r300');

keyboard
if false
    dimord={'delay corr','early corr', 'late corr', 'delay err', 'early err', 'late err'};
    save('ctd_accuracy.mat','dec_accu','c_reg','e_reg','dimord');
end



