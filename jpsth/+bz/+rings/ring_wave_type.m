function [out,sel]=ring_wave_type(in,opt)
arguments
    in
    opt.odor_only (1,1) logical
end
if opt.odor_only
    if all(ismember(in,[1 5]),'all') ...
            || all(ismember(in,[2 5]),'all') ...
            || all(ismember(in,[3 6]),'all') ...
            || all(ismember(in,[4 6]),'all')
        out='congru';
        sel='olf';

    elseif all(in==0,"all")
        out='nonmem';
        sel='none';
    elseif (any(in==5,"all") && any(in==6,"all")) ...
            ||(any(ismember(in,1:2),"all") && any(in==6,"all")) ...
            ||(any(ismember(in,3:4),"all") && any(in==5,"all")) ...
            ||(any(in==1 | in ==2,"all")+any(in==3 | in==4,"all")>1)
        out='incongru';
        sel='NA';
    else
        out='others';
        sel='NA';
    end
else
    if all(ismember(in,[1 5]),'all') ...
            || all(ismember(in,[2 5]),'all') ...
            || all(ismember(in,[3 6]),'all') ...
            || all(ismember(in,[4 6]),'all') ...
            || all(ismember(in,[1 7]),'all') ...
            || all(ismember(in,[2 8]),'all') ...
            || all(ismember(in,[3 7]),'all') ...
            || all(ismember(in,[4 8]),'all')
        out='congru';
        if any(ismember(in,1:4),'all')
            sel='both';
        elseif all(ismember(in,1:6),'all')
            sel='olf';
        elseif all(ismember(in,[1:4,7:8]),'all')
            sel='dur';
        else
            keyboard();
        end

    elseif all(in==0,"all")
        out='nonmem';
        sel='none';
    elseif (any(in==5,"all") && any(in==6,"all")) ...
            ||(any(in==7,"all") && any(in==8,"all")) ...
            ||(any(ismember(in,1:2),"all") && any(in==6,"all")) ...
            ||(any(ismember(in,3:4),"all") && any(in==5,"all")) ...
            ||(any(ismember(in,[1 3]),"all") && any(in==8,"all")) ...
            ||(any(ismember(in,[2 4]),"all") && any(in==7,"all")) ...
            ||(any(in==1,"all")+any(in==2,"all")+any(in==3,"all")+any(in==4,"all")>1)
        out='incongru';
        sel='NA';
    else
        out='others';
        sel='NA';
    end
end
end



