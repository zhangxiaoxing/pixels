function [is_diff,is_same]=diff_at_level(reg,level)
is_diff=false(numel(reg),1);
is_same=false(numel(reg),1);
for i=1:numel(reg)
    if ismember(reg{i}{1}{1},{'CH','BS'}) && ismember(reg{i}{2}{1},{'CH','BS'}) ...
            && ~isempty(reg{i}{1}{level}) && ~isempty(reg{i}{2}{level})
        if strcmp(reg{i}{1}{level},reg{i}{2}{level})
            is_same(i)=true;
        else
            is_diff(i)=true;
        end
    end
end
end
