keyboard()
denovo=false;
if denovo
    %wave_id
    meta=ephys.util.load_meta();
    waveids=ephys.get_wave_id(meta.sess,meta.allcid,early=true,ctx=false);
    %Dursel
    dur=wave.duration_wave('ctx',false,'extra_ITI',true);
    %Ordinal
    ord=wave.ordinal_MI(denovo=true,ctx=false);
    save('odor_dur_ordinal.mat','dur','meta','ord','waveids');
else
    load('odor_dur_ordinal.mat','dur','meta','ord','waveids');
end


odorsel=waveids>0;
bothsel=waveids>4;
dursel=any(dur.wrsp(:,1:7)<0.05,2);
ordisel=any(ord.anovanp(:,1:7)<0.05,2);

%ctx
ctxsel=strcmp(meta.reg_tree(2,:),'CTX').';

%nuclei
cnusel=strcmp(meta.reg_tree(2,:),'CNU').';

ctx_dist=typecount(ctxsel,odorsel,dursel,ordisel);
cnu_dist=typecount(cnusel,odorsel,dursel,ordisel);

ctx_prop=ctx_dist./sum(ctx_dist).*100;
cnu_prop=cnu_dist./sum(cnu_dist).*100;

fh=figure('Color','w');
bh=bar([ctx_prop;cnu_prop],'stacked');
legend(fliplr(bh),fliplr({'Non','Ordinal','Duration','Ord+Dur','Odor','Odor+Ord','Odor+Dur','Odor+Ord+Dur'}),'Location','eastoutside')
xlim([0.35,2.65])
set(gca,'XTick',1:2,'XTickLabel',{'Cortex','Nuclei'})
ylim([0,100])
ylabel('Fraction of neurons (%)');

both_dist=typecount(ctxsel&bothsel,odorsel,dursel,ordisel);
either_dist=typecount(ctxsel&odorsel&~bothsel,odorsel,dursel,ordisel);

both_prop=both_dist./sum(both_dist).*100;
either_prop=either_dist./sum(either_dist).*100;

fh=figure('Color','w');
bh=bar([both_prop;either_prop],'stacked');
legend(fliplr(bh),fliplr({'Non','Ordinal','Duration','Ord+Dur','Odor','Odor+Ord','Odor+Dur','Odor+Ord+Dur'}),'Location','eastoutside')
xlim([0.35,2.65])
set(gca,'XTick',1:2,'XTickLabel',{'Both 3s and 6s','Either 3s or 6s'})
ylim([0,100])
ylabel('Fraction of neurons');

function typecount=typecount(sel0,odorsel,dursel,ordisel)
typecount=nan(1,8);
for comb=0:7
    if bitand(comb,4)
        sel1=odorsel;
    else
        sel1=~odorsel;
    end
    if bitand(comb,2)
        sel2=dursel;
    else
        sel2=~dursel;
    end
    if bitand(comb,1)
        sel3=ordisel;
    else
        sel3=~ordisel;
    end
    typecount(comb+1)=nnz(sel0&sel1&sel2&sel3);
end
end