function sigtext(xx,yy,pp,bonf,color)
if ~exist('bonf','var')
    bonf=1;
end

if ~exist('color','var')
    color='k';
end

if pp<(0.001/bonf)
    text(xx,yy,'***','HorizontalAlignment','center','Color',color)
elseif pp<(0.01/bonf)
    text(xx,yy,'**','HorizontalAlignment','center','Color',color)
elseif pp<(0.05/bonf)
    text(xx,yy,'*','HorizontalAlignment','center','Color',color)
else
    text(xx,yy,'ns','HorizontalAlignment','center','Color',color)
end