% 3onlyS1:1,3onlyS2:2,6onlyS1:3,6onlyS2:4,bothS1:5,bothS2:6,nonmem:0,switch,etc:-1,
% require early and 6_late_half:-2
function out=get_wave_id(meta,opt)
arguments
    meta
    opt.split_6 (1,1) logical = false
end

assert(isa(meta,"struct"),'Function updated to receive meta structure as input')

sess=meta.sess;
cid=meta.allcid;

if isscalar(sess)
    sess=ones(size(cid)).*sess;
end
if numel(sess)~=numel(cid)
    throw(MException('getwave:arguments','Session and id number don''t match'))
end
persistent wave_id_map opt_
if isempty(wave_id_map) || ~isequaln(opt,opt_)
%     meta=ephys.util.load_meta();

    wave_id_map=containers.Map('KeyType','uint64','ValueType','int8');
    % 0, none; 1,3sec only;2,6sec only;3, both;-1,other
    for ii=1:numel(meta.sess)
        key=uint64(meta.sess(ii))*100000+uint64(meta.allcid(ii));

        if meta.mem_type_3(ii)==0 && meta.mem_type_6(ii)==0
            wave_id_map(key)=0;
        elseif ismember(meta.mem_type_3(ii),1:2) && meta.mem_type_6(ii)==0
            wave_id_map(key)=1;
        elseif ismember(meta.mem_type_3(ii),3:4) && meta.mem_type_6(ii)==0
            wave_id_map(key)=2;
        elseif ismember(meta.mem_type_6(ii),1:2) && meta.mem_type_3(ii)==0
            wave_id_map(key)=3;
        elseif ismember(meta.mem_type_6(ii),3:4) && meta.mem_type_3(ii)==0
            wave_id_map(key)=4;
        elseif ismember(meta.mem_type_3(ii),1:2) && ismember(meta.mem_type_6(ii),1:2)
            wave_id_map(key)=5;
        elseif ismember(meta.mem_type_3(ii),3:4) && ismember(meta.mem_type_6(ii),3:4)
            wave_id_map(key)=6;
        else
            wave_id_map(key)=-1;
        end
    end
end
key=uint64(sess).*100000+uint64(cid);
out=arrayfun(@(x) wave_id_map(x),key);
opt_=opt;
% out=cell2mat(wave_id_map.values(num2cell(key)));
end