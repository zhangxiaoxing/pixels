keyboard();
[sig,pair]=bz.load_sig_pair('pair',true);
% fh=figure('Color','w','Position',[32,32,400,225]);
meta=ephys.util.load_meta();
waveid=ephys.get_wave_id(meta.sess,meta.allcid);
anovameta=wave.get_dur_waveid();
%sel for wave id 7&8
waveid(waveid==0 & anovameta.dur_waveid==3)=7;
waveid(waveid==0 & anovameta.dur_waveid==6)=8;

for usess=reshape(anovameta.sess,1,[])
    %wave id 7&8
    asel=anovameta.sess==usess;
    ssel=sig.sess==usess;
    [~,loc]=ismember(sig.suid(ssel,:),anovameta.allcid(asel));
    sig.wave_id(ssel,:)=subsref(waveid(asel),struct(type='()',subs={{loc}}));

    psel=pair.sess==usess;
    [~,loc]=ismember(pair.suid(psel,:),anovameta.allcid(asel));
    pair.wave_id(psel,:)=subsref(waveid(asel),struct(type='()',subs={{loc}}));
end

[is_diff,is_same,~,~]=bz.util.diff_at_level(sig.reg,'hierarchy',false);
sigrsel=[is_diff(:,5),is_same(:,5),true(size(sig.sess))];

[is_diff,is_same,~,~]=bz.util.diff_at_level(pair.reg,'hierarchy',false);
pairrsel=[is_diff(:,5),is_same(:,5),true(size(pair.sess))];


ttl={'Different region','Same region','Overall'};
for rselidx=1:3
    hatmat=nan(9,9);
    cimat=nan(9,9,2);
    cntmat=nan(9,9,2);

    for lead=0:8
        for follow=0:8
            sigcnt=nnz(sigrsel(:,rselidx) & sig.wave_id(:,1)==lead & sig.wave_id(:,2)==follow);
            paircnt=nnz(pairrsel(:,rselidx) & pair.wave_id(:,1)==lead & pair.wave_id(:,2)==follow);
            [phat,pci]=binofit(sigcnt,paircnt);
            hatmat(lead+1,follow+1)=phat;
            cimat(lead+1,follow+1,:)=pci;
            cntmat(lead+1,follow+1,:)=[sigcnt,paircnt];
        end
    end
    hatnorm=hatmat-nnz(sigrsel(:,rselidx))./nnz(pairrsel(:,rselidx));
    fh=figure('Color','w');
    hold on
    imagesc(hatnorm,max(abs(hatnorm(:))).*[-1,1])
    plot([0.5,9.5],[0.5,9.5],'--w')
    colormap('turbo')
    colorbar();
    lbls={'Others','S1D3','S2D3','S1D6','S2D6','S1D3D6','S2D3D6','S1S2D3','S1S2D6'};
    set(gca(),'XTick',1:9,'XTickLabel',lbls,'YTick',1:9,'YTickLabel',lbls)
    xlabel('Lead');
    ylabel('Follow');
    title(ttl{rselidx})
%     truesize(gcf(),[360,360])
    exportgraphics(fh,'FC_assembly.pdf','ContentType','vector','Append',true);
end