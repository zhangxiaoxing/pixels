function out=fc_tag(in, antiCausal)
persistent warned
if isempty(warned)
    warned=false;
end
if ~exist('antiCausal','var')
    antiCausal=true;
    if ~warned
        disp('Will try anti-causal direction by default');
        warned=true;
        keyboard();
    end
    
end

out=in;
out(:,3:5)=0;

for i=1:size(in,1)-1
    if in(i,2)==1
        j=i+1;
        while in(j,2)==2 && j<size(in,1)
            if in(j,1)<=in(i,1)+0.01 && in(j,1)>in(i,1)+0.002
                out(i,3)=1;
                out(j,3)=1;
            end
            j=j+1;
        end
    end
end
if antiCausal
    for i=1:size(in,1)-1
        if in(i,2)==2
            j=i+1;
            while in(j,2)==1 && j<size(in,1)
                if in(j,1)<=in(i,1)+0.01 && in(j,1)>in(i,1)+0.002
                    out(i,4)=1;
                    out(j,4)=1;
                end
                j=j+1;
            end
        end
    end
end
for i=1:size(in,1)-1
    if in(i,2)==1
        j=i+1;
        while in(j,2)==2 && j<size(in,1)
            if in(j,1)<=in(i,1)+0.06 && in(j,1)>in(i,1)+0.052
                out(i,5)=1;
                out(j,5)=1;
            end
            j=j+1;
        end
    end
end
end