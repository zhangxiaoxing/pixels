% [bhat,bci]=binofit(both(1),both(2));keyboard();
[sig,pair]=bz.load_sig_pair('pair',true);
% fh=figure('Color','w','Position',[32,32,400,225]);
meta=ephys.util.load_meta();
waveid=ephys.get_wave_id(meta.sess,meta.allcid);
anovameta=ephys.selectivity_anova();
[dur_sense_mix,dur_exclu,~,sens_exclu]=ephys.get_dul_sel();


for usess=reshape(anovameta.sess,1,[])
    asel=anovameta.sess==usess;
    ssel=sig.sess==usess;
    [~,loc]=ismember(sig.suid(ssel,:),anovameta.allcid(asel));
    sig.wave_id(ssel,:)=subsref(waveid(asel),struct(type='()',subs={{loc}}));

    psel=pair.sess==usess;
    [~,loc]=ismember(pair.suid(psel,:),anovameta.allcid(asel));
    pair.wave_id(psel,:)=subsref(waveid(asel),struct(type='()',subs={{loc}}));
end

[is_diff,is_same,~,~]=bz.util.diff_at_level(sig.reg,'hierarchy',false);
sigrsel=[is_diff(:,5),is_same(:,5),is_diff(:,5)|is_same(:,5)];

[is_diff,is_same,~,~]=bz.util.diff_at_level(pair.reg,'hierarchy',false);
pairrsel=[is_diff(:,5),is_same(:,5),is_diff(:,5)|is_same(:,5)];

ttl={'Different region','Same region','Overall'};
%%
for rselidx=1:3
    both=[nnz(sigrsel(:,rselidx) & (all(sig.wave_id==5,2) | all(sig.wave_id==6,2))),...
        nnz(pairrsel(:,rselidx) & (all(pair.wave_id==5,2) | all(pair.wave_id==6,2)))];
    durOppoSens=[nnz(sigrsel(:,rselidx) & ((all(ismember(sig.wave_id,[1 2 5 6 7]),2) | all(ismember(sig.wave_id,[3 4 5 6 8]),2)) & ~(all(ismember(sig.wave_id,[1 3 5]),2) | all(ismember(sig.wave_id,[2 4 6]),2)))),...
        nnz(pairrsel(:,rselidx) &((all(ismember(pair.wave_id,[1 2 5 6 7]),2) | all(ismember(pair.wave_id,[3 4 5 6 8]),2)) & ~(all(ismember(pair.wave_id,[1 3 5]),2) | all(ismember(pair.wave_id,[2 4 6]),2))))];
    sensOppoDur=[nnz(sigrsel(:,rselidx) & ((all(ismember(sig.wave_id,[1 3 5]),2) | all(ismember(sig.wave_id,[2 4 6]),2)) & ~(all(ismember(sig.wave_id,[1 2 5 6 7]),2) | all(ismember(sig.wave_id,[3 4 5 6 8]),2)))),...
        nnz(pairrsel(:,rselidx) &((all(ismember(pair.wave_id,[1 3 5]),2) | all(ismember(pair.wave_id,[2 4 6]),2)) & ~(all(ismember(pair.wave_id,[1 2 5 6 7]),2) | all(ismember(pair.wave_id,[3 4 5 6 8]),2))))];
    oppoSensOppoDur=[nnz(sigrsel(:,rselidx) & (~(all(ismember(sig.wave_id,[1 2 5 6 7]),2) | all(ismember(sig.wave_id,[3 4 5 6 8]),2) | all(ismember(sig.wave_id,[1 3 5]),2) | all(ismember(sig.wave_id,[2 4 6]),2) | any(sig.wave_id==0,2)))),...
        nnz(pairrsel(:,rselidx) &(~(all(ismember(pair.wave_id,[1 2 5 6 7]),2) | all(ismember(pair.wave_id,[3 4 5 6 8]),2) | all(ismember(pair.wave_id,[1 3 5]),2) | all(ismember(pair.wave_id,[2 4 6]),2) | any(pair.wave_id==0,2))))];
    nonmem=[nnz(sigrsel(:,rselidx) & (all(sig.wave_id==0,2))),...
        nnz(pairrsel(:,rselidx) & (all(pair.wave_id==0,2)))];

    [bhat,bci]=binofit(both(1),both(2));
    [dhat,dci]=binofit(durOppoSens(1),durOppoSens(2));
    [shat,sci]=binofit(sensOppoDur(1),sensOppoDur(2));
    [xhat,xci]=binofit(oppoSensOppoDur(1),oppoSensOppoDur(2));
    [nhat,nci]=binofit(nonmem(1),nonmem(2));

    fh=figure('Color','w');
    hold on
    bh=bar([bhat,shat,dhat,nhat],'FaceColor','w');
    errorbar(bh.XEndPoints,bh.YEndPoints,...
        subsref([bci;sci;dci;nci],struct(type='()',subs={{':',1}})).'-bh.YEndPoints,...
        subsref([bci;sci;dci;nci],struct(type='()',subs={{':',2}})).'-bh.YEndPoints,...
        'k.');

    text(1,mean(ylim()),sprintf('%d / %d',both(1),both(2)),'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(3, max(ylim()),sprintf('%d / %d',durOppoSens(1),durOppoSens(2)),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle')
    text(2,max(ylim()),sprintf('%d / %d',sensOppoDur(1),sensOppoDur(2)),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle')
    %     text(4,mean(ylim()),sprintf('%d / %d',oppoSensOppoDur(1),oppoSensOppoDur(2)),'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(4,mean(ylim()),sprintf('%d / %d',nonmem(1),nonmem(2)),'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle')

    title(ttl{rselidx});
    set(gca(),'XTick',1:4,'XTickLabel',{'Same sensory-Same duration','Same sensory-Opposite duration','Same duration-Opposite sensory','Nonmemory'})

    [~,~,p]=crosstab([1*ones(both(2),1);2*ones(sensOppoDur(2),1);3*ones(durOppoSens(2),1);4*ones(nonmem(2),1)],...
        [(1:both(2))>both(1),(1:sensOppoDur(2))>sensOppoDur(1),(1:durOppoSens(2))>durOppoSens(1),(1:nonmem(2))>nonmem(1)]);
    [~,~,p1]=crosstab([1*ones(both(2),1);4*ones(nonmem(2),1)],...
        [(1:both(2))>both(1),(1:nonmem(2))>nonmem(1)]);
    [~,~,p2]=crosstab([2*ones(sensOppoDur(2),1);4*ones(nonmem(2),1)],...
        [(1:sensOppoDur(2))>sensOppoDur(1),(1:nonmem(2))>nonmem(1)]);
    [~,~,p3]=crosstab([3*ones(durOppoSens(2),1);4*ones(nonmem(2),1)],...
        [(1:durOppoSens(2))>durOppoSens(1),(1:nonmem(2))>nonmem(1)]);

    disp(ttl{rselidx})
    disp([p,p1,p2,p3])
    ylabel('Functional coupling rate');
end