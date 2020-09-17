fstr=cell(1,6);
for bin=1:6
        fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
% load('reg_keep.mat','reg_set')

rings4=cell(114,6);
parfor I=1:114
    %& diff(fstr{1}.reg_chain_S1,1,2)
    lbound=100000*I;
    ubound=100000*(I+1);
    for bin=1:6
        sel=fstr{bin}.conn_chain_S1(:,1)>lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
        if nnz(sel)>0
            onering=count_quadruplet(fstr{bin}.conn_chain_S1(sel,:),fstr{bin}.reg_chain_S1(sel,:),fstr{bin}.pref_chain_S1(sel,:),bin);
            if ~isempty(onering)
                onering(:,5)=bin;
                onering(:,6)=1;%sample
                rings4{I,bin}=onering;
            end
        end
    end
end
keyboard
save('rings.mat','rings4')
rings=cell(114,6);
parfor I=1:114
    %& diff(fstr{1}.reg_chain_S1,1,2)
    lbound=100000*I;
    ubound=100000*(I+1);
    for bin=1:6
        sel=fstr{bin}.conn_chain_S1(:,1)>lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
        if nnz(sel)>0
            rings{I,bin}=count_triplet(fstr{bin}.conn_chain_S1(sel,:),fstr{bin}.reg_chain_S1(sel,:),fstr{bin}.pref_chain_S1(sel,:),bin);
        end
    end
end
save('rings.mat','rings')


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

function out=count_quadruplet(in,reg,pref,bin)
out=[];
pre_unit_set=unique(in(:,1));
for i=reshape(pre_unit_set,1,[])
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
    for j=reshape(mono_post,1,[])
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
        for k=reshape(sec_post,1,[])
            third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
            ring=in(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,1);
            if ~isempty(ring)
                subgrp=repmat(ring,1,4);
                subgrp(:,1)=i;
                subgrp(:,2)=j;
                subgrp(:,3)=k;
                out=[out;subgrp];
            end
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



