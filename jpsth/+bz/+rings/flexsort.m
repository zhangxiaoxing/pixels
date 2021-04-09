function out=flexsort(in)
arguments
    in double %list of rings
end
out=in;
for i=1:size(out,1)
    [~,sess]=min(out(i,:));
    while sess>1
        out(i,:)=circshift(out(i,:),1);
        [~,sess]=min(out(i,:));
    end
end
end
