rpt=100;

for bin=1:6
    fbin(bin)=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
%    for bin=5
%        if ~isequal(fbin(bin).pair_chain,fbin(bin+1).pair_chain)
%            keyboard
%        end
%    end
load('reg_keep.mat');
% load('io_sel.mat');
greymatter=find(cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set));
% any(cell2mat(cellfun(@(x) x(:,1)<100,ioselstats,'UniformOutput',false)))
% nansel=io_entire_delay(:,1)<100;
% sel=find(~nansel & greymatter);
su_set=unique(fbin(1).conn_chain_S1(:,1));
input_to_count=nan(length(reg_set),6);
output_from_count=nan(length(reg_set),6);
recip_count=nan(length(reg_set),6);

input_to_count_shuf=nan(length(reg_set),6,100);
output_from_count_shuf=nan(length(reg_set),6,100);
recip_count_shuf=nan(length(reg_set),6,100);

parfor reg_target=1:length(reg_set)
    if ismember(reg_target,greymatter)
        disp(reg_target)
        for bin=1:6
            corrlist=[];
            for su=reshape(su_set,1,[])
                su2reg=nnz(fbin(bin).conn_chain_S1(:,1)==su & fbin(bin).reg_chain_S1(:,2)==reg_target);
                reg2su=nnz(fbin(bin).conn_chain_S1(:,2)==su & fbin(bin).reg_chain_S1(:,1)==reg_target);
                corrlist=[corrlist;su2reg,reg2su];
            end
            recip_count(reg_target,bin)=nnz(corrlist(:,1)>0 & corrlist(:,2)>0);
            output_from_count(reg_target,bin)=nnz(corrlist(:,2)>0);
            input_to_count(reg_target,bin)=nnz(corrlist(:,1)>0);
            
            for rr=1:rpt
                shuf=shuffle_conn_chain(fbin(bin).conn_chain_S1);
                corrlist=[];
                for su=reshape(su_set,1,[])
                    su2reg=nnz(shuf(:,1)==su & shuf(:,2)==reg_target);
                    reg2su=nnz(shuf(:,2)==su & shuf(:,1)==reg_target);
                    corrlist=[corrlist;su2reg,reg2su];
                end
                recip_count_shuf(reg_target,bin,rr)=nnz(corrlist(:,1)>0 & corrlist(:,2)>0);
                output_from_count_shuf(reg_target,bin,rr)=nnz(corrlist(:,2)>0);
                input_to_count_shuf(reg_target,bin,rr)=nnz(corrlist(:,1)>0);
            end

        end
        %keyboard
    end
end
%keyboard
%save('su_region_recip.mat','per_reg_recip_idx') 
save('su_region_recip.mat','recip_count','input_to_count','output_from_count','recip_count_shuf','input_to_count_shuf','output_from_count_shuf')
keyboard
load('su_region_recip.mat')
wing=diff(ioselstats{1}(:,[3 7]),1,2);
per_reg_recip_idx=[per_reg_recip_idx,(1:140)',wing];
per_reg_recip_idx=per_reg_recip_idx(~isnan(per_reg_recip_idx(:,1)),:);

[~,idces]=sort(per_reg_recip_idx(:,1));
fh=figure('Color','w','Position',[0,100,1600,400]);
bar(1:length(idces),per_reg_recip_idx(idces,1))
set(gca,'XTick',1:length(idces),'XTickLabel',reg_set(per_reg_recip_idx(idces,2)),'XTickLabelRotation',90)
ylim([0,1])
ylabel('reciprocal fraction');
if ~verLessThan('matlab','9.8')
    exportgraphics(fh,'su_region_recip.pdf')
end

figure('Color','w','Position',[100,100,400,300]);
scatter(per_reg_recip_idx(:,1),per_reg_recip_idx(:,3),'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
xlabel('reciprocal fraction');
ylabel('WING');
[r,p]=corr(per_reg_recip_idx(:,1),per_reg_recip_idx(:,3));
legend(sprintf('r=%.3f,p=%.3f',r,p));
if ~verLessThan('matlab','9.8')
    exportgraphics(fh,'su_region_recip_WING.pdf')
end

function out=shuffle_conn_chain(in)
out=in;
for i=1:114
    lbound=i*100000;
    ubound=(i+1)*100000;
    sel=in(:,1)>lbound & in(:,1)<ubound;
    if nnz(sel)>0
        seldata=in(sel,:);
        shufsel=randperm(nnz(sel));
        shufdata=seldata(shufsel,2);
        out(sel,2)=shufdata;        
    end
end
end
