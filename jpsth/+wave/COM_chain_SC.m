function COM_chain_SC(opt)
arguments
    opt.peak (1,1) logical = false
    opt.strict (1,1) logical = true %strict FC criteria
end
load('sums_conn.mat','sums_conn_str');
warning('partial iteration for illustration')
for fidx=87%1:numel(sums_conn_str)
    disp(fidx);
    if opt.strict
        if isempty(sums_conn_str(fidx).ccg_sc)
            continue
        end
        ccg=sums_conn_str(fidx).ccg_sc;
        strict_sel=ccg(:,3)>0 & ccg(:,6)<20 & ccg(:,5)<20 & ccg(:,4)<350 & ccg(:,4)>252;
        oneccg=ccg(strict_sel,:);
    else
        oneccg=sums_conn_str(fidx).ccg_sc;
    end
    
    onecom=wave.get_com_map('onepath',sums_conn_str(fidx).folder,'peak',opt.peak);
    skey=fieldnames(onecom);
    for samp=["s1","s2"]
        mapkeys=cell2mat(onecom.(skey{1}).(samp).keys);
        typesel=all(ismember(int32(oneccg(:,1:2)),mapkeys),2);
        typesigcon=oneccg(typesel,1:2);
        con_com_diff=arrayfun(@(x) onecom.(skey{1}).(samp)(x),int32(oneccg(typesel,1:2)));
        % disp(min(con_com_diff,[],'all'))
        % disp(max(con_com_diff,[],'all'))
        if opt.peak
            dirsel=con_com_diff(:,2)>=con_com_diff(:,1);
        else
            dirsel=con_com_diff(:,2)-con_com_diff(:,1)>(50/250);
        end
        %         dircom=con_com_diff(dirsel,:);
        dirsigcon=typesigcon(dirsel,:);
        upre=unique(dirsigcon(:,1)).';
        
        chains=cell(0);
        for i=upre
            onechain=cell(0);
            cpre=i;
            while true
                newpair=dirsigcon(ismember(dirsigcon(:,1),cpre),:);
                if isempty(newpair)
                    if numel(onechain)>2
                        chains{end+1}=onechain;
                    end
                    break
                else
                    onechain{end+1}={...
                        newpair,...
                        arrayfun(@(x) onecom.(skey{1}).(samp)(x),newpair)};
                    cpre=newpair(:,2);
                end
            end
        end
        if numel(chains)>0
            for ii=1:numel(chains)
                onechain=chains{ii};
                if ~(dt(onechain)>4)
                    continue
                end
                fh=figure('Color','w','Position',[100,100,250,250]);
                hold on;
                yidx=1;
                text(onechain{1}{2}(1,1),yidx,num2str(onechain{1}{1}(1,1)));
                for jj=1:numel(onechain)
                    plot(onechain{jj}{2}.',repmat([yidx;yidx+1],[1,size(onechain{jj}{2},1)]),'ko','MarkerFaceColor','k');
                    arrayfun(@(x) ...
                        quiver(onechain{jj}{2}(x,1),yidx,diff(onechain{jj}{2}(x,:)),1,'-r',...
                        'MaxHeadSize',2./sqrt(diff(onechain{jj}{2}(x,:)).^2+1)...
                        ),1:size(onechain{jj}{2},1));
                    arrayfun(@(x) ...
                        text(onechain{jj}{2}(x,2),yidx+1,num2str(onechain{jj}{1}(x,2))),...
                        1:size(onechain{jj}{2},1));                    
                    yidx=yidx+1;
                    
                end
                ylim([0.5,yidx+0.5]);
                xlim([0,24]);
                set(gca(),'XTick',[0,8,16,24],'XTickLabel',[0,2,4,6],'YTick',0:numel(onechain)+1)
                xlabel('COM time in delay (s)')
                ylabel('Functional coupling order in wave ->');
                title(sprintf('Session %s, sample %s',skey{1},samp));
                fhccg=figure('Color','w','Position',[100,100,250,500]);
                vlen=3%size(onechain,2);
                hlen=2%max(cellfun(@(x) size(x{1},1), onechain));
                cidx=1;
                for kk=1:vlen
                    for ll=1:size(onechain{kk}{1},1)
                        conn=onechain{kk}{1}(ll,:);
%                         disp(conn)
                        subplot(vlen,hlen,cidx);cidx=cidx+1;%(kk-1)*hlen+ll);
                        ccgsel=oneccg(:,1)==conn(1) & oneccg(:,2)==conn(2);
                        plot(oneccg(ccgsel,7:end),'-r');
                        arrayfun(@(x) xline(x,'--k'),[226,251,276]);
                        xlim([191,311]);
                        set(gca(),'XTick',201:25:301,'XTickLabel',-20:10:20)
                        text(191,max(ylim()),num2str(conn(1)),'HorizontalAlignment','left','VerticalAlignment','top');
                        text(301,max(ylim()),num2str(conn(2)),'HorizontalAlignment','right','VerticalAlignment','top');
                    end
                end
                keyboard();
                a
            end
        else
            disp('Empty Chain');
        end
    end
end
end
function out=dt(onechain)
out=max(cellfun(@(x) max(x{2},[],'all'),onechain))...
    - min(cellfun(@(x) min(x{2},[],'all'),onechain));
end