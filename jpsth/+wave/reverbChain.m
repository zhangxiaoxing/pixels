sig_str=my.load_sig_pair();

onesample(sig_str.s1,'prefsamp',1,'plot',true);
onesample(sig_str.s1,'prefsamp',2,'plot',true);

onesample(sig_str.s2,'prefsamp',2,'plot',true);
onesample(sig_str.s2,'prefsamp',1,'plot',true);
% onesample(sig);


function onesample(sig,opt)
arguments
    sig (1,1) struct
    opt.plot (1,1) logical = false
    opt.prefsamp (1,1) int32 = 1
end

    prefsamp=opt.prefsamp;

    bin1postsel=(sig.per_bin(:,1,1)==prefsamp & all(sig.per_bin(:,1:2,2)==prefsamp,2));
    
    bin2presel=post2pre(sig,bin1postsel);
    
    bin2postsel=bin2presel & all(sig.per_bin(:,2:3,2)==prefsamp,2);
    
    bin3presel=post2pre(sig,bin2postsel);
    
    bin3postsel=bin3presel & all(sig.per_bin(:,3:4,2)==prefsamp,2);
    
    bin4presel=post2pre(sig,bin3postsel);

    bin4postsel=bin4presel & all(sig.per_bin(:,4:5,2)==prefsamp,2);
    
    bin5presel=post2pre(sig,bin4postsel);
    
    bin5postsel= bin5presel & all(sig.per_bin(:,5:6,2)==prefsamp,2);
    
    
%%
if opt.plot
    cmap=[1,1,1;1,0,0];
    fh=figure('Color','w','Position',[32,32,200,9000]);
    subplot(6,1,1)
    imagesc(sig.per_bin(bin1postsel,1:6,1))
    colormap(cmap);
    ax1=gca();
    subplot(6,1,2)
    imagesc(sig.per_bin(bin1postsel,1:6,2))
    colormap(cmap);
    ax2=gca();
    subplot(6,1,3)
    imagesc(sig.per_bin(bin2postsel,1:6,2))
    colormap(cmap);
    ax3=gca();
    subplot(6,1,4)
    imagesc(sig.per_bin(bin3postsel,1:6,2))
    colormap(cmap);
    ax4=gca();
    subplot(6,1,5)
    imagesc(sig.per_bin(bin4postsel,1:6,2))
    colormap(cmap);
    ax5=gca();
    subplot(6,1,6)
    imagesc(sig.per_bin(bin5postsel,1:6,2))
    colormap(cmap);
    ax6=gca();

    for ax=[ax1 ax2 ax3 ax4 ax5 ax6]
        ax.XTick=[];
        ax.YTick=[1,max(ax.YTick)];
    end
    ax6.XTick=1:6;
    ax6.XLabel.String='Time bin (sec)';
%     keyboard
    exportgraphics(fh,sprintf('prolonged_showcase_%d.pdf',prefsamp),'ContentType','vector');
end
end

function out=post2pre(sig,sel)
    out=false(size(sel));
    id=sig.suid(sel,2);
    sess=sig.sess(sel);
    usess=reshape(unique(sess),1,[]);
    for s=usess
        sid=id(sess==s);
        out(ismember(sig.suid(:,1),sid) & sig.sess==s)=true;
    end
end
