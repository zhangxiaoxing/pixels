function out=getTCOM(sess,cid)
arguments
    sess (:,1) double
    cid (:,1) double
end

if numel(sess)~=numel(cid)
    throw(MException('get_coactive_profile:arguments','Session and id number don''t match'))
end
persistent TCOM_map
if isempty(TCOM_map)
    com_map=wave.get_com_map('curve',false,'wave','any','delay',6);
    TCOM_map=containers.Map('KeyType','uint64','ValueType','double');

    % 0, none; 1,3sec only;2,6sec only;3, both;-1,other
    for fn=reshape(fieldnames(com_map),1,[])
        fs=fn{1};
        sessnum=str2double(replace(fs,'s',''));
        key_s1=num2cell(uint64(sessnum)*100000+uint64(cell2mat(com_map.(fs).s1.keys)));
        if isempty(key_s1), continue;end
        raw_s1=cell2mat(com_map.(fs).s1.values);
        sessmap_s1=containers.Map(key_s1,raw_s1);
        TCOM_map=[TCOM_map;sessmap_s1];
        key_s2=num2cell(uint64(sessnum)*100000+uint64(cell2mat(com_map.(fs).s2.keys)));
        if isempty(key_s2), continue;end
        raw_s2=cell2mat(com_map.(fs).s2.values);
        sessmap_s2=containers.Map(key_s2,raw_s2);
        TCOM_map=[TCOM_map;sessmap_s2];
    end
end
keys=uint64(sess).*100000+uint64(cid);
keysel=TCOM_map.isKey(num2cell(keys));
out=nan(size(sess));
out(keysel)=cell2mat(TCOM_map.values(num2cell(keys(keysel))));
end


function debug
    ephys.getTCOM(sig.sess,sig.suid(:,1))
end