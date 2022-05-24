function reverbChainBZ(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.prefix (1,:) char = 'BZWT'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
sig=bz.load_sig_pair('type',opt.type,'prefix',opt.prefix);
onesample(sig,'prefsamp',2,'plot',true);
% onesample(sig);
end

function onesample(sig,opt)
arguments
    sig (1,1) struct
    opt.plot (1,1) logical = false
    opt.prefsamp (1,1) int32 = 1
end
 
    prefsamp=opt.prefsamp;

    bin1sel1=(sig.per_bin(:,1,1)==prefsamp & all(sig.per_bin(:,1:2,2)==prefsamp,2));
    
    bin2presel=post2pre(sig,bin1sel1);
    
    bin2sel1=bin2presel & all(sig.per_bin(:,2:3,2)==prefsamp,2);
    
    bin3presel=post2pre(sig,bin2sel1);
    
    bin3sel1=bin3presel & all(sig.per_bin(:,3:4,2)==prefsamp,2);
    
    bin4presel=post2pre(sig,bin3sel1);

    bin4sel1=bin4presel & all(sig.per_bin(:,4:5,2)==prefsamp,2);
    
    bin5presel=post2pre(sig,bin4sel1);
    
    bin5sel1= bin5presel & all(sig.per_bin(:,5:6,2)==prefsamp,2);
    
    
%%
if opt.plot
    cmap=[1,1,1;1,0,0];
    fh=figure('Color','w','Position',[100,100,200,600]);
    subplot(6,1,1)
    imagesc(sig.per_bin(bin1sel1,1:6,1))
    colormap(cmap);
    ax1=gca();
    subplot(6,1,2)
    imagesc(sig.per_bin(bin1sel1,1:6,2))
    colormap(cmap);
    ax2=gca();
    subplot(6,1,3)
    imagesc(sig.per_bin(bin2sel1,1:6,2))
    colormap(cmap);
    ax3=gca();
    subplot(6,1,4)
    imagesc(sig.per_bin(bin3sel1,1:6,2))
    colormap(cmap);
    ax4=gca();
    subplot(6,1,5)
    imagesc(sig.per_bin(bin4sel1,1:6,2))
    colormap(cmap);
    ax5=gca();
    subplot(6,1,6)
    imagesc(sig.per_bin(bin5sel1,1:6,2))
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

% conn_pref_mat=[fstr{1}.pref_chain_S1(bin1sel,1:6);...
%     nan(25,6);...
%     fstr{1}.pref_chain_S1(bin1sel,7:12);...
%     nan(25,6);...
%     fstr{2}.pref_chain_S1(bin2sel,7:12);...
%     nan(25,6);...
%     fstr{3}.pref_chain_S1(bin3sel,7:12);...
%     nan(25,6);...
%     fstr{4}.pref_chain_S1(bin4sel,7:12);...
%     nan(25,6);...
%     fstr{5}.pref_chain_S1(bin5sel,7:12)];
% % close all
% fh=figure();
% imagesc(conn_pref_mat);
% % exportgraphics(fh,'enhance_showcase_v01.pdf','ContentType','vector')
