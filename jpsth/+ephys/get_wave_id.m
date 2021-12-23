% 3onlyS1:1,3onlyS2:2,6onlyS1:3,6onlyS2:4,bothS1:5,bothS2:6,nonmem:0,switch,etc:-1,
% require early and 6_late_half:-2
function out=get_wave_id(sess,cid,opt)
arguments
    sess double
    cid double
    opt.early (1,1) logical = true
end
if isscalar(sess)
    sess=ones(size(cid)).*sess;
end
if numel(sess)~=numel(cid)
    throw(MException('getwave:arguments','Session and id number don''t match'))
end
persistent wave_id_map early_
if isempty(wave_id_map) || opt.early~=early_
    meta3=ephys.util.load_meta('delay',3);
    meta6=ephys.util.load_meta('delay',6);
    % all(meta3.sess==meta6.sess,"all")
    % all(meta3.allcid==meta6.allcid,"all")
    wave_id_map=containers.Map('KeyType','uint64','ValueType','int8');
    % 0, none; 1,3sec only;2,6sec only;3, both;-1,other
    for ii=1:numel(meta3.sess)
        key=uint64(meta3.sess(ii))*100000+uint64(meta3.allcid(ii));
        if meta3.mem_type(ii)==0 && meta6.mem_type(ii)==0
            wave_id_map(key)=0;        
        elseif ismember(meta3.mem_type(ii),1:2) && meta6.mem_type(ii)==0
            wave_id_map(key)=1;
        elseif ismember(meta3.mem_type(ii),3:4) && meta6.mem_type(ii)==0
            wave_id_map(key)=2;
        elseif ismember(meta6.mem_type(ii),1:2) && meta3.mem_type(ii)==0
            if ~opt.early || (opt.early && any(meta6.per_bin(1:3,ii)==1))
                wave_id_map(key)=3;
            else
                wave_id_map(key)=-2;
            end
        elseif ismember(meta6.mem_type(ii),3:4) && meta3.mem_type(ii)==0
            if ~opt.early || (opt.early && any(meta6.per_bin(1:3,ii)==2))
                wave_id_map(key)=4;
            else
                wave_id_map(key)=-2;
            end
        elseif ismember(meta3.mem_type(ii),1:2) && ismember(meta6.mem_type(ii),1:2)
            if ~opt.early || (opt.early && any(meta6.per_bin(1:3,ii)==1))
                wave_id_map(key)=5;
            else
                wave_id_map(key)=-2;
            end
        elseif ismember(meta3.mem_type(ii),3:4) && ismember(meta6.mem_type(ii),3:4)
            if ~opt.early || (opt.early && any(meta6.per_bin(1:3,ii)==2))
                wave_id_map(key)=6;
            else
                wave_id_map(key)=-2;
            end
        else
            wave_id_map(key)=-1;
        end
    end
end
key=uint64(sess).*100000+uint64(cid);
out=arrayfun(@(x) wave_id_map(x),key);
early_=opt.early;
% out=cell2mat(wave_id_map.values(num2cell(key)));
end