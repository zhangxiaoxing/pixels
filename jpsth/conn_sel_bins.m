findIdx=false;
if findIdx
    subtotal=sum(cat(3,conn_mat_all{:}),3);
    [~,idx]=max(subtotal(:));
    [row,col]=ind2sub(size(subtotal),idx);
end
% col=23,xaxis,source ,'DP', row=77,yaxis,target,'TTd', 

for bin=1:6
    disp(bin);
    load(sprintf('XCORR_stats_delay_6_%d_%d_2msbin.mat',bin,bin+1));
    for scell=stats
        s=scell{1};
%         if ((strcmp(s.reg_su1,'AId') && strcmp(s.reg_su2,'AON') && s.s1_peak_significant && s.AIs1<0) ||... %AIs <0 means 1->2
%             (strcmp(s.reg_su1,'AON') && strcmp(s.reg_su2,'AId') && s.s1_peak_significant && s.AIs1>0)) &&... already using selective neurons?
%             (any(s.prefered_sample_su1(2:end)) && any(s.prefered_sample_su2(2:end))) %selective doesn't mean favour
        if s.s1_peak_significant && any(s.prefered_sample_su1(2:end)) && any(s.prefered_sample_su2(2:end))
            key=s.fileidx*1000*1000+s.su1_label_idx*1000+s.su2_label_idx;
            idx=find(su_var_t(:,1)==key);
            if s.prefered_sample_su1(bin+1)>0 && s.prefered_sample_su2(bin+1)>0
                selType=2;
            else
                selType=1;
            end
            if isempty(idx)
                su_var_t(end+1,:)=zeros(1,7);
                su_var_t(end,1)=key;
                su_var_t(end,bin+1)=selType;
            else
                su_var_t(idx,bin+1)=selType;
            end
        end
    end
end


for bin=2:6
    d_conn_pair(bin)=nnz((su_var_t(:,bin)>0) ~= (su_var_t(:,bin+1)>0))/length(su_var_t);
    d_conn_num(bin)=(nnz(su_var_t(:,bin+1)>0)-nnz(su_var_t(:,bin)>0))/length(su_var_t);
    
    d_sel_pair(bin)=nnz((su_var_t(:,bin)>1) ~= (su_var_t(:,bin+1)>1))/length(su_var_t);
    d_sel_num(bin)=(nnz(su_var_t(:,bin+1)>1)-nnz(su_var_t(:,bin)>1))/length(su_var_t);
end


figure('Color','w','Position',[20,20,600,800])
subplot(1,2,1);
imagesc(su_var_t(:,2:end));
cmap=[1 1 1;0 0 1;1 0 0];
colormap(cmap)
xlabel('delay bin (s)')
ylabel('SU #')
% title('DP -> TTd')

subplot(3,2,2);
hold on
ph1=plot(sum(su_var_t(:,2:end)>0),'b-','LineWidth',1.5);
ph2=plot(sum(su_var_t(:,2:end)>1),'r-','LineWidth',1.5);
ylim([0,max(ylim())])
xlim([1,6])
xlabel('delay time bin (s)')
ylabel('connection pairs')
legend([ph1,ph2],{'connective','selective'});


subplot(3,2,4)
hold on
phn=plot(2:6,d_conn_num(2:end),'-b','LineWidth',1.5);
php=plot(2:6,d_conn_pair(2:end),'--b','LineWidth',1.5);
legend([phn,php],{'delta global connectivity','delta pair connectivity'})
xlabel('delay bin (s)')
ylabel('fraction of all pairs')


subplot(3,2,6)
hold on
phn=plot(2:6,d_sel_num(2:end),'-r','LineWidth',1.5);
php=plot(2:6,d_sel_pair(2:end),'--r','LineWidth',1.5);
legend([phn,php],{'delta global selectivity','delta pair selectivity'})
xlabel('delay bin (s)')
ylabel('fraction of all pairs')

sgtitle('brain-wide')
print('conn_dyn_brain-wide.png','-dpng')

% colorbar()


key=su_var_t(:,1);
source=idivide(int32(key),int32(100));
target=rem(key,100);
stable=all(su_var_t(:,2:end)');

stable_sel=all(su_var_t(:,2:end)'>1);


unisource=unique(source);

stable_count=zeros(0,2);
transient_count=zeros(0,2);
for s=unisource'
    stable_count(end+1,:)=[s,nnz(source==s & stable')];
    transient_count(end+1,:)=[s,nnz(source==s & ~stable')];
end


[shist,edges]=histcounts(stable_count(:,2),1:2:60);
thist=histcounts(transient_count(:,2),1:2:60);
figure('Color','w','Position',[100,100,450,300])
bh=bar(2:2:59, [thist',shist'],'stacked');
legend({'transient','sustained'})
xlabel('number of target neuron from same source neuron')
ylabel('number of source neuron')
print('Neuron-neuron-hist.png','-dpng')


figure()
histogram(stable_count(:,2),1:2:80)
figure()
histogram(transient_count(:,2),1:2:80)


load('join_reg_set.mat')

