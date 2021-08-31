[sig,pair]=bz.load_sig_pair();
rpts=50;
sums=cell(116,2);
parfor sess=1:2
    ffcc=[];
    nonfc=[];
    disp(sess)
    sess_fc=sig.suid(sig.sess==sess,:);
    fc_set=[bitshift(sess_fc(:,1),14)+sess_fc(:,2);bitshift(sess_fc(:,2),14)+sess_fc(:,1)];
    for samp=1:2
        [stats,key_map]=shuffle_one(sess,rpts,sprintf('s%da',samp),sprintf('s%db',samp));
        comvar=cellfun(@(x) nanstd(x),stats);
        for key=reshape(cell2mat(key_map.keys),1,[])
            suid1=double(bitshift(key,-14));
            suid2=double(bitand(key,hex2dec('3FFF')));
            if ismember(key,fc_set)
                ffcc=[ffcc;sess,samp,suid1,suid2,comvar(key_map(key))];
            else
                nonfc=[nonfc;sess,samp,suid1,suid2,comvar(key_map(key))];
            end
        end
    end
    sums(sess,:)={ffcc,nonfc};
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
    one_comb=nchoosek(cell2mat(samp_key),2);
    for ci=1:size(one_comb,1)
        if ~key_map.isKey(getkey(one_comb(ci,:)))
            key_map(getkey(one_comb(ci,:)))=statsIdx;
            stats{statsIdx}=nan(1,100);
            statsIdx=statsIdx+1;
        end
        stats{key_map(getkey(one_comb(ci,:)))}(rpt)=com_map.(fs).(field_a)(one_comb(ci,1))-com_map.(fs).(field_a)(one_comb(ci,2));
        stats{key_map(getkey(one_comb(ci,:)))}(rpt+1)=com_map.(fs).(field_b)(one_comb(ci,1))-com_map.(fs).(field_b)(one_comb(ci,2));
    end
end
end
