%featii,epochii,jj,double(glmxmeta(comb2(jj,1),:)),double(glmxmeta(comb2(jj,2),:)),mdl.Coefficients.Estimate.',mdl.Rsquared.Ordinary,mdl.ModelCriterion.AICc
keyboard()
t=two_reg_corr_list;
ureg=unique([t(:,6);t(:,9)]);
minrdist=[];
for io1=1:2
    for io2=1:2
        iosel=t(:,4)==io1 & t(:,7)==io2;
        for fii=1:6
            fsel=t(:,1)==fii;
            for uii=1:size(ureg,1)
                if rem(uii,10)==0
                    disp(uii)
                end
                u1sel=t(:,6)==ureg(uii);
                u2reg=unique(t(fsel & iosel & u1sel,9));
                for r2=reshape(u2reg,1,[])
                    rowsel=iosel & fsel & u1sel & t(:,9)==r2;
                    minr=t(rowsel,14).';
                    if numel(minr)~=6
                        keyboard()
                    end
                    minrdist=[minrdist;fii,io1,io2,ureg(uii),r2,minr,min(minr)];
                end
            end
        end
    end
end
%%
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
postfixmap=containers.Map({1,2},{'From ','To '});
feat_l={'Mixed-sample','Exclusive-sample','Mixed-duration','Exclusive-duration','Mixed-position','Exclusive-position'};
mm=min(minrdist(:,5:end),[],2);
pid=find(mm>0.36);
%%
for pid=231481%reshape(pid,1,[])
e1id=find(t(:,1)==minrdist(pid,1) & t(:,2)==1 & t(:,4)==minrdist(pid,2) & t(:,7)==minrdist(pid,3) & t(:,6)==minrdist(pid,4)  & t(:,14)==minrdist(pid,5));
r2=t(e1id,9);

cross_epoch=find(t(:,1)==minrdist(pid,1) & t(:,4)==minrdist(pid,2) & t(:,7)==minrdist(pid,3) & t(:,6)==minrdist(pid,4) & t(:,9)==r2);
coeffs=t(cross_epoch,[2,10:15]);
fh=figure('Color','w','Position',[32,32,1280,400]);
subplot(1,2,1)
plot(coeffs(:,1),sqrt(coeffs(:,6)))
xlim([0.5,6.5])
ylim([0,1])
set(gca(),'XTick',1:6,'XTickLabel',{'Sample','EarlyDelay','LateDelay','Test','PostReward','PreNextTrial'})
ylabel('Pearson''s r')
title(feat_l{minrdist(pid,1)});

subplot(1,2,2)
hold on
ih=plot(coeffs(:,1),coeffs(:,2),'--k');
r1h=plot(coeffs(:,1),coeffs(:,3),'-r');
r2h=plot(coeffs(:,1),coeffs(:,4),'-b');
inth=plot(coeffs(:,1),coeffs(:,5),'-m');
xlim([0.5,6.5])
set(gca(),'XTick',1:6,'XTickLabel',{'Sample','EarlyDelay','LateDelay','Test','PostReward','PreNextTrial'})
ylabel('Coefficient')
r1r=idmap.ccfid2reg(minrdist(pid,4));
r2r=idmap.ccfid2reg(r2);
r1p=postfixmap(minrdist(pid,2));
r2p=postfixmap(minrdist(pid,3));
title(sprintf('GLM %d',pid));
legend([r1h,r2h,inth,ih],{[r1p,r1r{1}],[r2p,r2r{1}],'Interaction','Intercept'},'Location','north','Orientation','horizontal')
exportgraphics(fh,sprintf('two_reg_glm_%d.pdf',pid),'ContentType','vector')
end