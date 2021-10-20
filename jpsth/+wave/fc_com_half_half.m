[sig,pair]=bz.load_sig_pair();
rpts=50;
sums=cell(116,3);
parfor sess=1:116
    [fc_fwd,fc_rev,non_fc]=deal([]);
    disp(sess)
    sess_fc=sig.suid(sig.sess==sess,:);
    com_map=wave.get_com_map('onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',false,'rnd_half',false);
    fc_set_fwd=bitshift(sess_fc(:,1),14)+sess_fc(:,2);
    fc_set_rev=bitshift(sess_fc(:,2),14)+sess_fc(:,1);
    for samp=1:2
        [stats,key_map]=shuffle_one(sess,rpts,sprintf('s%da',samp),sprintf('s%db',samp));
        comvar=cellfun(@(x) nanstd(x),stats);
        commm=cellfun(@(x) nanmean(x),stats);
        for key=reshape(cell2mat(key_map.keys),1,[])
            suid1=double(bitshift(key,-14));
            suid2=double(bitand(key,hex2dec('3FFF')));
            if ~all(arrayfun(@(x) com_map.(sprintf('s%d',sess)).(sprintf('s%d',samp)).isKey(x),[suid1,suid2]))
                continue;
            end
            fc_com=arrayfun(@(x) com_map.(sprintf('s%d',sess)).(sprintf('s%d',samp))(x),[suid1,suid2]);
            if (ismember(key,fc_set_fwd) && diff(fc_com)>0) || (ismember(key,fc_set_rev) && diff(fc_com)<0)
                fc_fwd=[fc_fwd;sess,samp,suid1,suid2,commm(key_map(key)),comvar(key_map(key))];
            elseif (ismember(key,fc_set_rev) && diff(fc_com)>0) || (ismember(key,fc_set_fwd) && diff(fc_com)<0)
                fc_rev=[fc_rev;sess,samp,suid1,suid2,commm(key_map(key)),comvar(key_map(key))];
            else
                non_fc=[non_fc;sess,samp,suid1,suid2,commm(key_map(key)),comvar(key_map(key))];
            end
        end
    end
    sums(sess,:)={fc_fwd,fc_rev,non_fc};
end
save('fc_com_half_half.mat','sums')

function plot()
load('fc_com_half_half.mat','sums');
fc_fwd=cell2mat(sums(:,1));
fc_rev=cell2mat(sums(:,2));
non_fc=cell2mat(sums(:,3));
fc_fwd(:,7)=(fc_fwd(:,6).*0.25).^2;
fc_rev(:,7)=(fc_rev(:,6).*0.25).^2;
non_fc(:,7)=(non_fc(:,6).*0.25).^2;
mm=[mean(fc_fwd(:,7)),mean(fc_rev(:,7)),mean(non_fc(:,7))];
ci_fwd=bootci(500,@(x) mean(x),fc_fwd(:,7));
ci_rev=bootci(500,@(x) mean(x),fc_rev(:,7));
boots=bootstrp(500,@(x) mean(x),non_fc(:,7));
ci_non=[prctile(boots,2.5);prctile(boots,97.5)];
fh=figure('Color','w','Position',[32,32,190,190]);
hold on
bar(mm,'FaceColor','none','EdgeColor','k','LineWidth',1);
errorbar(1:3,mm,[ci_fwd(1),ci_rev(1),ci_non(1)]-mm,[ci_fwd(2),ci_rev(2),ci_non(2)]-mm,'k.','CapSize',15)
ylabel('TCOM latency variance(sec2)')
set(gca(),'XTick',1:3,'XTickLabel',{'Prog. F.C.','Regres. F.C.','No coupling'},'XTickLabelRotation',30);
xlim([0.4,3.6]);
exportgraphics(fh,'FC_com_half_half.pdf');
end


function [stats,key_map]=shuffle_one(sid,rpts,field_a,field_b)
getkey=@(x) bitshift(x(1),14)+x(2);
stats=cell(0);
statsIdx=1;
key_map=containers.Map('KeyType','int32','ValueType','int32');
fs=sprintf('s%d',sid);
for rpt=1:2:2*rpts
    com_map=wave.get_com_map('onepath',['SPKINFO/',ephys.sessid2path(sid)],'curve',false,'rnd_half',true);
    samp_key=num2cell(intersect(cell2mat(com_map.(fs).(field_a).keys),cell2mat(com_map.(fs).(field_b).keys)));
    if numel(samp_key)<2
        break;
    end
    one_comb=nchoosek(cell2mat(samp_key),2);
    for ci=1:size(one_comb,1)
        if ~key_map.isKey(getkey(one_comb(ci,:)))
            key_map(getkey(one_comb(ci,:)))=statsIdx;
            stats{statsIdx}=nan(1,100);
            statsIdx=statsIdx+1;
        end
        stats{key_map(getkey(one_comb(ci,:)))}(rpt)=com_map.(fs).(field_a)(one_comb(ci,2))-com_map.(fs).(field_a)(one_comb(ci,1));
        stats{key_map(getkey(one_comb(ci,:)))}(rpt+1)=com_map.(fs).(field_b)(one_comb(ci,2))-com_map.(fs).(field_b)(one_comb(ci,1));
    end
end
end
