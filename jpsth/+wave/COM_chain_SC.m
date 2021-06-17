function COM_chain_SC(opt)
arguments
    opt.peak (1,1) logical = false
    opt.strict (1,1) logical = true %strict FC criteria
end
load('sums_conn.mat','sums_conn_str');
warning('partial iteration for illustration')
for fidx=1:numel(sums_conn_str)
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

        if opt.peak
            dirsel=con_com_diff(:,2)>con_com_diff(:,1)...
                &con_com_diff(:,2)<=con_com_diff(:,1)+8;
        else
            dirsel=con_com_diff(:,2)-con_com_diff(:,1)>(50/250);
        end
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
                split_chains=wave.recursive_chain(chains{ii},[]);
                %TODO split chains with shared components
                for pp=1:2:size(split_chains,1)
                    onechain=split_chains(pp:pp+1,:);
                    if onechain(2,end)-onechain(2,1)<4
                        continue
                    end
                    fh=figure('Color','w','Position',[100,100,500,500]);
                    subplot(2,2,3);
                    hold on;
                    for jj=1:size(onechain,2)
                        text(onechain(2,jj),jj,num2str(onechain(1,jj)));
                        plot(onechain(2,jj),jj,'ko','MarkerFaceColor','k');
                    end
                    for jj=1:size(onechain,2)-1
                        plot(onechain(2,jj:jj+1),[jj;jj+1],'r-');
                    end
                    ylim([0.5,size(onechain,2)+0.5]);
                    xlim([0,24]);
                    set(gca(),'XTick',[0,8,16,24],'XTickLabel',[0,2,4,6],'YTick',0:numel(onechain)+1)
                    xlabel('COM time in delay (s)')
                    ylabel('Functional coupling order in wave ->');
                    title(sprintf('Session %s, sample %s',skey{1},samp));
                    
                    
                    for kk=1:size(onechain,2)-1
                        conn=onechain(1,kk:kk+1);
                        subplot(2,size(onechain,2)-1,kk);
                        ccgsel=oneccg(:,1)==conn(1) & oneccg(:,2)==conn(2);
                        plot(oneccg(ccgsel,7:end),'-r');
                        arrayfun(@(x) xline(x,'--k'),[226,251,276]);
                        xlim([191,311]);
                        set(gca(),'XTick',201:25:301,'XTickLabel',-20:10:20)
                        text(191,max(ylim()),num2str(conn(1)),'HorizontalAlignment','left','VerticalAlignment','top');
                        text(301,max(ylim()),num2str(conn(2)),'HorizontalAlignment','right','VerticalAlignment','top');
                    end
                    
                    
                    subplot(2,2,4)
                    hold on
                    suids=onechain(1,:);
%                     pdim=ceil(sqrt(numel(suids)));
                    cmap=colormap('cool');
                    deltac=floor(size(cmap,1)./numel(suids));
                    for mm=1:numel(suids)
    %                     subplot(pdim,pdim,mm);

                        plot(smooth(onecom.(skey{1}).(sprintf('%scurve',samp))(suids(mm)),5),...
                            '-','Color',cmap((mm-1)*deltac+1,:));
    %                     title(suids(mm));
                        xline(onecom.(skey{1}).(samp)(suids(mm)),'--','Color',cmap((mm-1)*deltac+1,:));
                    end
                    ylim([0,0.3])
                    set(gca(),'XTick',0:8:24,'XTickLabel',0:2:6)
                    keyboard();
                end
            end
        else
            disp('Empty Chain');
        end
    end
end
end
