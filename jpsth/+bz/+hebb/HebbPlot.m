% based on miyashita's trial-shuffle pipeline
% Not in active use. 2023.09.25

%sample=1;
load('rings.mat','rings'); %from ring_list.m dimord=(msize, session, bin, sample)
hebbPattern=cell(0,4);
for bin=1:6
    fstr{bin}=load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
for b3=1:6
    r3all=cell2mat(rings(1,:,b3,sample)');
    for b4=b3:min(b3+1,6)
        r4all=cell2mat(rings(2,:,b4,sample)');
        r4all=r4all(:,1:4);
        for r3idx=1:length(r3all)
            for r4idx=1:length(r4all)
                if nnz(ismember(r4all(r4idx,:),r3all(r3idx,:)))==2
                    r4place=find(ismember(r4all(r4idx,:),r3all(r3idx,:)));
                    if isequal(r4place,[1 3])...
                            || isequal(r4place,[2 4])...
                        r3place=find(ismember(r3all(r3idx,:),r4all(r4idx,:)));
                        switch sum(r3place)    
                            case 4 %1 3 or 3 1
                                final_ring124=r3all(r3idx,1:3);
%                                 post15=
                            case 3 % 1 2
                                final_ring124=r3all(r3idx,[2 3 1]);
                            case 5 % 2 3
                                final_ring124=r3all(r3idx,[3 1 2]);
                        end
                        for b4p=b4:min(b4+1,6)
                            r4pall=cell2mat(rings(2,:,b4p,sample)');
                            r4pall=r4pall(:,1:4);
                            for b4pidx=1:length(r4pall)
                                if nnz(ismember(r4pall(b4pidx,:),final_ring124))==3 ...
                                        && isequal(flexsort(r4pall(b4pidx,ismember(r4pall(b4pidx,:),final_ring124))),flexsort(final_ring124))
                                    pre14Idx=find(r4pall(b4pidx,:)==final_ring124(3));
                                    if (pre14Idx==1 && ~ismember(r4pall(b4pidx,4),final_ring124)) ...
                                            || (pre14Idx~=1 && ~ismember(r4pall(b4pidx,pre14Idx-1),final_ring124))
                                        post15idx=find(r4all(r4idx,:)==final_ring124(1))+1;
                                        if post15idx>4
                                            post15idx=post15idx-4;
                                        end
                                        if any(fstr{b4p}.conn_chain_S1(:,1)==final_ring124(1) & fstr{b4p}.conn_chain_S1(:,2)==r4all(r4idx,post15idx))
%                                             disp(r3all(r3idx,:))
%                                             disp(r4all(r4idx,:))
%                                             disp(r4pall(b4pidx,:))
                                            hebbPattern(end+1,:)={r3all(r3idx,:),r4all(r4idx,:),r4pall(b4pidx,:),[b3,b4,b4p]};
%                                             keyboard
                                        end
                                    end
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if sample==1
    hebbPatternS1=hebbPattern;
    save('hebb_pattern.mat','hebbPatternS1','-append');
elseif sample==2
    hebbPatternS2=hebbPattern;
    save('hebb_pattern.mat','hebbPatternS2','-append');
end

return
keyboard
fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
load('reg_keep.mat','reg_set')

inter_reg_cong_count=zeros(114,1);
sess_reg=cell(114,1);
for i=1:114
    lbound=100000*i;
    ubound=100000*(i+1);
    sel=fstr{6}.conn_chain_S1(:,1)>lbound & fstr{6}.conn_chain_S1(:,1)<ubound;
    inter_reg_cong_count(i)=nnz(diff(fstr{6}.reg_chain_S1(sel,:),1,2) & fstr{6}.pref_chain_S1(sel,6)==fstr{6}.pref_chain_S1(sel,12));
    sess_reg{i}=arrayfun(@(x) reg_set{x}, unique(fstr{6}.reg_chain_S1(sel,:)),'UniformOutput',false);
end
[C,I]=sort(inter_reg_cong_count,'descend');
sess_reg=[C,sess_reg(I)];




%%
sel1=find(diff(fstr{1}.reg_chain_S1,1,2)...
    & diff(fstr{1}.pref_chain_S1(:,[1 7]),1,2)==0 ...
    & fstr{1}.pref_chain_S1(:,1)>0 ...
    & diff(fstr{1}.pref_chain_S1(:,6:7),1,2)==0);
for b1=reshape(sel1,[],1)
    post1=fstr{1}.conn_chain_S1(b1,2);
    sel2=find( ...
        fstr{2}.conn_chain_S1(:,1)==post1 ...
        & diff(fstr{2}.reg_chain_S1,1,2)...
        & diff(fstr{2}.pref_chain_S1(:,[1 7]),1,2)==0 ...
        & fstr{2}.pref_chain_S1(:,1)>0 ...
        & diff(fstr{2}.pref_chain_S1(:,6:7),1,2)==0);
    if isempty(sel2)
        continue
    end
%     for post2=reshape
end





%% combining rings, not working as expected
rings=cell(114,6);
parfor I=1:114
%& diff(fstr{1}.reg_chain_S1,1,2)
lbound=100000*I;
ubound=100000*(I+1);

for bin=1:6
    sel=fstr{bin}.conn_chain_S1(:,1)>lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
    rings{I,bin}=count_triplet(fstr{bin}.conn_chain_S1(sel,:),fstr{bin}.reg_chain_S1(sel,:),fstr{bin}.pref_chain_S1(sel,:),bin);
end
end


nnz(ismember(rings{4}(:),rings{5}(:)))
cc=fstr{5}.conn_chain_S1(sel,:);
nnz(ismember(cc(:,1),rings{4}(:)) & ismember(cc(:,2),rings{5}(:))) %indirect connn from 4 to 5


% unique(rings{3}(ismember(rings{3}(:),rings{4}(:))))

hebbnet=[];
for b4=1:size(rings{4},1)
    for b3=1:size(rings{3},1)
        if ~any(ismember(rings{3}(b3,:),rings{4}(b4,:)))
            continue
        end
        for b2=1:size(rings{2},1)
            if ~any(ismember(rings{2}(b2,:),rings{3}(b3,:)))
                continue
            end
            for b1=1:size(rings{1},1)
                if ~any(ismember(rings{1}(b1,:),rings{2}(b2,:)))
                    continue
                end
%                 keyboard
                hebbnet=[hebbnet;rem([rings{1}(b1,:),rings{2}(b2,:),rings{3}(b3,:),rings{4}(b4,:)],100000)];
            end
        end
    end
end

nodesCount=arrayfun(@(x) numel(unique(hebbnet(x,:))),1:length(hebbnet));
[nc,ni]=sort(nodesCount);
hebbnet=hebbnet(ni,:);
save('hebbplot.mat','rings','hebbnet','nc')






prefsamp=1;
geplist=[];
reglist=[];
bin1sel=(fstr{1}.pref_chain_S1(:,1)==prefsamp & all(fstr{1}.pref_chain_S1(:,7:8)==prefsamp,2));
% bin1postu=unique(fstr{1}.conn_chain_S1(bin1sel,2));
geplist=[geplist;fstr{1}.conn_chain_S1(bin1sel,:),ones(nnz(bin1sel),1),fstr{1}.peaks1(bin1sel)'];
reglist=[reglist;reshape(fstr{1}.conn_chain_S1(bin1sel,:),[],1),reshape(fstr{1}.reg_chain_S1(bin1sel,:),[],1)];


function out=count_triplet(in,reg,pref,bin)
out=[];
pre_unit_set=unique(in(:,1));
for i=reshape(pre_unit_set,1,[])
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
    for j=reshape(mono_post,1,[])
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
        ring=in(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,1);
        if ~isempty(ring)
            subgrp=[ring,ring,ring];
            subgrp(:,1)=i;
            subgrp(:,2)=j;
            out=[out;subgrp];
        end
    end
end
out=unique(flexsort(out),'rows');
end

function out=flexsort(in)
    out=in;
    for i=1:size(out,1)
        [~,I]=min(out(i,:));
        while I>1
            out(i,:)=circshift(out(i,:),1);
            [~,I]=min(out(i,:));
        end
    end
end

