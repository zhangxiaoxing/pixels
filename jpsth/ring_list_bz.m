%Caution: will not check for existance due to parfor limitation

%TODO check WM region before finalize

if (~exist('delay_data','var')) || (~exist('delay_inact_data','var')) ...
        ||(~exist('delay_shuf','var')) || (~exist('delay_shuf_inact','var'))
    disp('Missing necessary switch')
    return
end

fstr=cell(1,6);
bzthres=250;
for bin=1:6
    fstr{bin}=load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    s1sel=all(fstr{bin}.bz_spk_count_S1>bzthres,2);
    s2sel=all(fstr{bin}.bz_spk_count_S2>bzthres,2);
    fstr{bin}.bz_conn_reg_S1=fstr{bin}.bz_conn_reg_S1(s1sel,:);
    fstr{bin}.bz_conn_reg_S2=fstr{bin}.bz_conn_reg_S2(s2sel,:);
    fstr{bin}.bz_conn_chain_S1=fstr{bin}.bz_conn_chain_S1(s1sel,:);
    fstr{bin}.bz_conn_chain_S2=fstr{bin}.bz_conn_chain_S2(s2sel,:);
    fstr{bin}.bz_pref_S1= fstr{bin}.bz_pref_S1(s1sel,:);
    fstr{bin}.bz_pref_S2= fstr{bin}.bz_pref_S2(s2sel,:);
end
bin=-2;
fbase=load(sprintf('0906_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));


% load('reg_keep.mat','reg_set')
if (exist('delay_data','var') && delay_data) || (exist('delay_inact_data','var') && delay_inact_data)
    rings=cell(3,114,6,2);
    rings_inact=cell(3,114,6,2);
    for I=1:114
        for midx=1:3
            msize=midx+2;
            disp(I)
            lbound=100000*I;
            ubound=100000*(I+1);
            for bin=1:6
                sel11=fstr{bin}.bz_conn_chain_S1(:,1)>=lbound & fstr{bin}.bz_conn_chain_S1(:,1)<ubound & diff(fstr{bin}.bz_conn_reg_S1,1,2);
                if nnz(sel11)>0
                    if  delay_data
                        onering1=count_motif_congru(fstr{bin}.bz_conn_chain_S1(sel11,:),fstr{bin}.bz_conn_reg_S1(sel11,:),fstr{bin}.bz_pref_S1(sel11,:),bin,msize,1);
                        onering1=unique(flexsort(onering1),'rows');
                        sel12=fstr{bin}.bz_conn_chain_S2(:,1)>=lbound & fstr{bin}.bz_conn_chain_S2(:,1)<ubound & diff(fstr{bin}.bz_conn_reg_S2,1,2);
                        onering2=count_motif_congru(fstr{bin}.bz_conn_chain_S2(sel12,:),fstr{bin}.bz_conn_reg_S2(sel12,:),fstr{bin}.bz_pref_S2(sel12,:),bin,msize,2);
                        onering2=unique(flexsort(onering2),'rows');
                        rings(midx,I,bin,:)={onering1,onering2};
                    end
                    if  delay_inact_data
                        onering1=count_motif_congru_inact(fstr{bin}.bz_conn_chain_S1(sel11,:),fstr{bin}.bz_conn_reg_S1(sel11,:),fstr{bin}.bz_pref_S1(sel11,:),msize,1);
                        onering1=unique(flexsort(onering1),'rows');
                        sel12=fstr{bin}.bz_conn_chain_S2(:,1)>=lbound & fstr{bin}.bz_conn_chain_S2(:,1)<ubound & diff(fstr{bin}.bz_conn_reg_S2,1,2);
                        onering2=count_motif_congru_inact(fstr{bin}.bz_conn_chain_S2(sel12,:),fstr{bin}.bz_conn_reg_S2(sel12,:),fstr{bin}.bz_pref_S2(sel12,:),msize,2);
                        onering2=unique(flexsort(onering2),'rows');
                        rings_inact(midx,I,bin,:)={onering1,onering2};
                    end
                end
            end
            
        end
    end
    if exist('delay_data','var') && delay_data
        save('rings_bz.mat','rings','-append');
    end
    if exist('delay_inact_data','var') && delay_inact_data
        save('rings_bz.mat','rings_inact','-append');
    end
end


if (exist('delay_shuf','var') && delay_shuf) || (exist('delay_shuf_inact','var') && delay_shuf_inact)
    shufrpt=100;
    rings_shuf=cell(shufrpt,3,114,6,2);
    rings_shuf_inact=cell(shufrpt,3,114,6,2);
    for rpt=1:shufrpt
        disp(rpt)
        for bin=1:6
            [shufchainS1,shufregS1,shufprefS1]=shuffle_conn_chain(fstr{bin}.bz_conn_chain_S1,fstr{bin}.pair_chain,fstr{bin}.pair_reg,fstr{bin}.pref_pair);
            [shufchainS2,shufregS2,shufprefS2]=shuffle_conn_chain(fstr{bin}.bz_conn_chain_S2,fstr{bin}.pair_chain,fstr{bin}.pair_reg,fstr{bin}.pref_pair);
            parfor I=1:114
                for midx=1:3
                    msize=midx+2;
                    lbound=100000*I;
                    ubound=100000*(I+1);
                    sel21=shufchainS1(:,1)>=lbound & shufchainS1(:,1)<ubound & diff(shufregS1,1,2);
                    if nnz(sel21)>0
                        if delay_shuf
                            onering1=count_motif_congru(shufchainS1(sel21,:),shufregS1(sel21,:),shufprefS1(sel21,:),bin,msize,1);
                            onering1=unique(flexsort(onering1),'rows');
                            sel22=shufchainS2(:,1)>=lbound & shufchainS2(:,1)<ubound & diff(shufregS2,1,2);
                            onering2=count_motif_congru(shufchainS2(sel22,:),shufregS2(sel22,:),shufprefS2(sel22,:),bin,msize,2);
                            onering2=unique(flexsort(onering1),'rows');
                            rings_shuf(rpt,midx,I,bin,:)={onering1,onering2};
                        end
                        if  delay_shuf_inact
%                             disp('Go')
                            onering1=count_motif_congru_inact(shufchainS1(sel21,:),shufregS1(sel21,:),shufprefS1(sel21,:),msize,1);
                            onering1=unique(flexsort(onering1),'rows');
                            sel22=shufchainS2(:,1)>=lbound & shufchainS2(:,1)<ubound & diff(shufregS2,1,2);
                            onering2=count_motif_congru_inact(shufchainS2(sel22,:),shufregS2(sel22,:),shufprefS2(sel22,:),msize,2);
                            onering2=unique(flexsort(onering1),'rows');
                            rings_shuf_inact(rpt,midx,I,bin,:)={onering1,onering2};
                        end                        
                    end
                end
            end
        end
    end
    if exist('delay_shuf','var') && delay_shuf
        save('rings_bz.mat','rings_shuf','-append');
    end
    if exist('delay_shuf_inact','var') && delay_shuf_inact
        save('rings_bz.mat','rings_shuf_inact','-append');
    end
end

%% baseline
if exist('base_data','var') && base_data
    base_rings=cell(3,114,2);
    parfor I=1:114
        for midx=1:3
            msize=midx+2;
            disp(I)
            lbound=100000*I;
            ubound=100000*(I+1);
            sel11=fbase.bz_conn_chain_S1(:,1)>=lbound & fbase.bz_conn_chain_S1(:,1)<ubound & diff(fbase.bz_conn_reg_S1,1,2);
            if nnz(sel11)>0
                onering1=count_motif_congru_inact(fbase.bz_conn_chain_S1(sel11,:),fbase.bz_conn_reg_S1(sel11,:),fbase.bz_pref_S1(sel11,:),msize,1);
                onering1=unique(flexsort(onering1),'rows');
                sel12=fbase.bz_conn_chain_S2(:,1)>=lbound & fbase.bz_conn_chain_S2(:,1)<ubound & diff(fbase.bz_conn_reg_S2,1,2);
                onering2=count_motif_congru_inact(fbase.bz_conn_chain_S2(sel12,:),fbase.bz_conn_reg_S2(sel12,:),fbase.bz_pref_S2(sel12,:),msize,2);
                onering2=unique(flexsort(onering2),'rows');
                base_rings(midx,I,:)={onering1,onering2};
            end
        end
    end
    save('rings_bz.mat','base_rings','-append');
end
if exist('base_shuf','var') && base_shuf
    shufrpt=1000;
    base_rings_shuf=cell(shufrpt,3,114,2);
    for rpt=1:shufrpt
        disp(rpt)
        [shufchainS1,shufregS1,shufprefS1]=shuffle_conn_chain(fbase.bz_conn_chain_S1,fbase.pair_chain,fbase.pair_reg,fbase.pref_pair);
        [shufchainS2,shufregS2,shufprefS2]=shuffle_conn_chain(fbase.bz_conn_chain_S2,fbase.pair_chain,fbase.pair_reg,fbase.pref_pair);
        parfor I=1:114
            for midx=1:3
                msize=midx+2;
                lbound=100000*I;
                ubound=100000*(I+1);
                sel21=shufchainS1(:,1)>=lbound & shufchainS1(:,1)<ubound & diff(shufregS1,1,2);
                if nnz(sel21)>0
                    onering1=count_motif_congru_inact(shufchainS1(sel21,:),shufregS1(sel21,:),shufprefS1(sel21,:),msize,1);
                    onering1=unique(flexsort(onering1),'rows');
                    sel22=shufchainS2(:,1)>=lbound & shufchainS2(:,1)<ubound & diff(shufregS2,1,2);
                    onering2=count_motif_congru_inact(shufchainS2(sel22,:),shufregS2(sel22,:),shufprefS2(sel22,:),msize,2);
                    onering2=unique(flexsort(onering1),'rows');
                    base_rings_shuf(rpt,midx,I,:)={onering1,onering2};
                end
            end
        end
    end
    save('rings_bz.mat','base_rings_shuf','-append');
end




return


function out=count_motif_congru(in,reg,pref,bin,msize,sample)
out=[];
pre_unit_set=unique(in(:,1));
for i=pre_unit_set(:)'
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
    for j=mono_post(:)'
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
        if msize==3
            ring=in(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([i j],size(subgrp,1),1);
                out=[out;subgrp];
            end
        elseif msize==4
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                ring=in(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                if ~isempty(ring)
                    subgrp=[ring,ring,ring,ring];
                    subgrp(:,1:3)=repmat([i j k],size(subgrp,1),1);
                    out=[out;subgrp];
                end
            end
        elseif msize==5
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                third_post(third_post==i)=[];
                for l=third_post(:)'
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                    fourth_post(ismember(fourth_post,[i j]))=[];
                    ring=in(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring,ring];
                        subgrp(:,1:4)=repmat([i j k l],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                end
            end
        end
    end
end


end

function out=count_motif_congru_inact(in,reg,pref,msize,sample)
out=[];
pre_unit_set=unique(in(:,1));
for i=pre_unit_set(:)'
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
    for j=mono_post(:)'
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
        if msize==3
            ring=in(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([i j],size(subgrp,1),1);
                out=[out;subgrp];
            end
        elseif msize==4
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                ring=in(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                if ~isempty(ring)
                    subgrp=[ring,ring,ring,ring];
                    subgrp(:,1:3)=repmat([i j k],size(subgrp,1),1);
                    out=[out;subgrp];
                end
            end
        elseif msize==5
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                third_post(third_post==i)=[];
                for l=third_post(:)'
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                    fourth_post(ismember(fourth_post,[i j]))=[];
                    ring=in(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring,ring];
                        subgrp(:,1:4)=repmat([i j k l],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                end
            end
        end
    end
end


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


function [outconn,outreg,outpref]=shuffle_conn_chain(in,inpair,inpair_reg,inpair_pref)
outconn=nan(size(in));
outreg=nan(size(in));
outpref=nan(size(in,1),12);
combined=false(size(in,1),1);
for i=1:114
    lbound=i*100000;
    ubound=(i+1)*100000;
    sel=in(:,1)>=lbound & in(:,1)<ubound;
    combined=combined | sel;
    if nnz(sel)>0
        selpair=find(inpair(:,1)>=lbound & inpair(:,1)<ubound);
        shufsel=randperm(nnz(selpair));
        %        keyboard
        shufdata=inpair(selpair(shufsel(1:nnz(sel))),:);
        flipsel=randi(2,size(shufdata,1),1)>1;
        shufdata(flipsel,:)=shufdata(flipsel,[2,1]);
        outconn(sel,:)=shufdata;
        if exist('inpair_reg','var')
            outreg(sel,:)=inpair_reg(selpair(shufsel(1:nnz(sel))),:);
        end
        
        if exist('inpair_pref','var')
            outpref(sel,:)=inpair_pref(selpair(shufsel(1:nnz(sel))),:);
        end
    end
end
% sum(combined)
end
