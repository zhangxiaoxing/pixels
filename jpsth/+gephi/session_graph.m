
function session_graph(opt)
arguments
    opt.gen_data (1,1) logical = false
    opt.plot_scatter (1,1) logical = false
    opt.plot_bar (1,1) logical = false
end
if gen_data
    [sig,~]=bz.load_sig_pair('pair',false);
    meta=ephys.util.load_meta();
    ctxsel=all(squeeze(sig.reg(:,2,:))==688,2);
    for sess=reshape(unique(sig.sess),1,[])
        
        %meta select data region size
        
        [s1c,s2c,nonc,anyc]=deal({'Source','Target'});
        [s1n,s2n,nonn,anyn]=deal({'Id','Label'});
        
        sess_sel=find(sig.sess==sess & ctxsel);
        %same s1, same s2, nonmem, any
        % {source,target}
        for si=reshape(sess_sel,1,[])
            anyc=[anyc;{sig.suid(si,1),sig.suid(si,2)}];
            if all(sig.mem_type(si,:)==0,'all') % non
                nonc=[nonc;{sig.suid(si,1),sig.suid(si,2)}];
            elseif all(ismember(sig.mem_type(si,:),1:2),'all')
                s1c=[s1c;{sig.suid(si,1),sig.suid(si,2)}];
            elseif all(ismember(sig.mem_type(si,:),3:4),'all')
                s2c=[s2c;{sig.suid(si,1),sig.suid(si,2)}];
            end
        end
        writecell(s1c,fullfile('gephidata',sprintf('s1_edge_%03d.csv',sess)));
        writecell(s2c,fullfile('gephidata',sprintf('s2_edge_%03d.csv',sess)));
        writecell(nonc,fullfile('gephidata',sprintf('non_edge_%03d.csv',sess)));
        writecell(anyc,fullfile('gephidata',sprintf('any_edge_%03d.csv',sess)));
        
        s1n=[s1n;num2cell(repmat(meta.allcid(meta.sess==sess & ismember(meta.mem_type,1:2).' & strcmp(meta.reg_tree(2,:),'CTX').'),1,2))];
        s2n=[s2n;num2cell(repmat(meta.allcid(meta.sess==sess & ismember(meta.mem_type,3:4).' & strcmp(meta.reg_tree(2,:),'CTX').'),1,2))];
        nonn=[nonn;num2cell(repmat(meta.allcid(meta.sess==sess & (meta.mem_type==0).' & strcmp(meta.reg_tree(2,:),'CTX').'),1,2))];
        anyn=[anyn;num2cell(repmat(meta.allcid(meta.sess==sess & (meta.mem_type>=0).' & strcmp(meta.reg_tree(2,:),'CTX').'),1,2))];
        
        writecell(s1n,fullfile('gephidata',sprintf('s1_node_%03d.csv',sess)));
        writecell(s2n,fullfile('gephidata',sprintf('s2_node_%03d.csv',sess)));
        writecell(nonn,fullfile('gephidata',sprintf('non_node_%03d.csv',sess)));
        writecell(anyn,fullfile('gephidata',sprintf('any_node_%03d.csv',sess)));
        %     keyboard()
    end
end
if opt.plot_scatter || opt.plot_bar
    congru_S1=readmatrix(fullfile('gephidata','s1_gephi_graph_sums.csv'));
    congru_S2=readmatrix(fullfile('gephidata','s2_gephi_graph_sums.csv'));
    nonmem=readmatrix(fullfile('gephidata','non_gephi_graph_sums.csv'));
    any=readmatrix(fullfile('gephidata','any_gephi_graph_sums.csv'));
    node_sel=congru_S1(:,1)>20 & nonmem(:,1)>20;
    if opt.plot_scatter
        plotOne(congru_S1(node_sel,[1 3]),nonmem(node_sel,[1 3]),'Connected comp.');
        plotOne(congru_S1(node_sel,[1 4]),nonmem(node_sel,[1 4]),'Avg. cluster coef.');
        plotOne(congru_S1(node_sel,[1 5]),nonmem(node_sel,[1 5]),'Density');
        plotOne(congru_S1(node_sel,[1 6]),nonmem(node_sel,[1 6]),'Avg. degree');
        plotOne(congru_S1(node_sel,[1 7]),nonmem(node_sel,[1 7]),'Avg. path length');
    end
    if opt.plot_bar
        plotOneBar(congru_S1(node_sel,3),nonmem(node_sel,3),'Connected comp.');
        plotOneBar(congru_S1(node_sel,4),nonmem(node_sel,4),'Avg. cluster coef.');
        plotOneBar(congru_S1(node_sel,5),nonmem(node_sel,5),'Density');
        plotOneBar(congru_S1(node_sel,6),nonmem(node_sel,6),'Avg. degree');
        plotOneBar(congru_S1(node_sel,7),nonmem(node_sel,7),'Avg. path length');
        
        plotOneBar(congru_S2(node_sel,3),nonmem(node_sel,3),'Connected comp.');
        plotOneBar(congru_S2(node_sel,4),nonmem(node_sel,4),'Avg. cluster coef.');
        plotOneBar(congru_S2(node_sel,5),nonmem(node_sel,5),'Density');
        plotOneBar(congru_S2(node_sel,6),nonmem(node_sel,6),'Avg. degree');
        plotOneBar(congru_S2(node_sel,7),nonmem(node_sel,7),'Avg. path length');
    end
end
end

function plotOne(exp,ref,lbl)
figure()
hold on;
scatter(exp(:,1),exp(:,2),4,'r');
scatter(ref(:,1),exp(:,2),4,'b');
title(lbl);
end


function fh=plotOneBar(congru,incongru,lbl,opt)
arguments
    congru
    incongru
    lbl (1,:) char
    opt.exportgraph (1,1) logical = false
end
fh=figure('Color','w','Position',[100,100,125,200]);
hold on;
ydata=[incongru,congru];
xdata=ones(size(ydata));
if contains(lbl,'comp.')
    xdata(:,1)=xdata(:,1)+rand(size(ydata,1),1)*0.3;
    xdata(:,2)=xdata(:,2)+1-rand(size(ydata,1),1)*0.3;
else
    xdata(:,1)=xdata(:,1)+rand(size(ydata,1),1)*0.2;
    xdata(:,2)=xdata(:,2)+1-rand(size(ydata,1),1)*0.2;
end
plot(xdata',ydata','-k')
plot(xdata(:,1),ydata(:,1),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','none')
plot(xdata(:,2),ydata(:,2),'ko','MarkerSize',4,'MarkerFaceColor','none','MarkerEdgeColor','k')
errorbar([0.75,2.25],mean(ydata),std(ydata)/sqrt(size(ydata,1)),'ko','MarkerSize',6);
xlim([0.5,2.5])
ylim([0,max(ylim())])
ylabel(lbl)
xlabel(sprintf('%0.4f',signrank(congru,incongru)));
set(gca,'XTick',[])
if contains(lbl,'comp.')
    set(gca,'YScale','log')
    ylim([0.9,11])
    ax=gca;
    ax.YAxis.MinorTickValues=1:9;
    ax.YAxis.MinorTick='on';
    ax.TickLength=[0.03,0.06];
end
if opt.exportgraph
    exportgraphics(fh,sprintf('func_conn_stats_%s.pdf',lbl),'ContentType','vector');
end
end


