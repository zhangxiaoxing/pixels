function fh=cmap_histogram(data,win,cmap)
fh=figure('Color','w');
hold on;
for ii=1:(numel(win)-1)
    ycount=nnz(data>win(ii) & data<=win(ii+1));
    fill([win(ii),win(ii),win(ii+1),win(ii+1)],...
        [0,ycount,ycount,0],...
        cmap(ii,:),'EdgeColor','k');

end
end

