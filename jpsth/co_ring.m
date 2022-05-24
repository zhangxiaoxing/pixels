sessidx=46;
if ~exist('rings','var')
    load rings.mat
end
ringall=[];
pref_all=[];
for midx=1:3
r1=[];
for bb=1:6
    r1=[r1;cell2mat(rings(midx,sessidx,bb,1)')];
end
r2=[];
for bb=1:6
    r2=[r2;cell2mat(rings(midx,sessidx,bb,2)')];
end
r1(:,midx+3)=1;r2(:,midx+3)=2;ringsm=[r1(:,1:(midx+2));r2(:,1:(midx+2))];samps=[r1(:,midx+3);r2(:,midx+3)];
if size(ringsm,2)<5
   ringsm(:,end+1:5)=-1; 
end
ringall=[ringall;ringsm];
pref_all=[pref_all;samps];
end


ringsu=unique(ringall,'rows');
comat=[];
for i=1:(length(rings1)-1)
    disp(i);
    for j=(i+1):length(rings1)
        r1=rings1(i,:);
        r2=rings1(j,:);
        r1=r1(r1>0);
        r2=r2(r2>0);
        comat(i,j)=nnz(ismember(r1,r2));
    end
end
comatBk=comat;
for i=2:(length(rings1)-1)
    for j=1:(i-1)
        comat(i,j)=comat(j,i);
    end
end

close all
fh=figure('Color','w','Position',[100,100,255,225]);
imagesc(comat,[0 4])
cb=colorbar();
cb.Label.String='Shared neuron';
xlabel('Ring #')
ylabel('Ring #')
exportgraphics(fh,'rings_shared_neurons.pdf')
