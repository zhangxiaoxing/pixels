function [is_diff,is_same]=diff_at_level(reg)
is_diff=false(numel(reg),6);
is_same=false(numel(reg),6);
for i=1:numel(reg)
    for treedep=1:6
        if ismember(reg{i}{1}{1},{'CH','BS'}) && ismember(reg{i}{2}{1},{'CH','BS'}) ...
                && ~isempty(reg{i}{1}{treedep}) && ~isempty(reg{i}{2}{treedep})
            if strcmp(reg{i}{1}{treedep},reg{i}{2}{treedep})
                is_same(i,treedep)=true;
            else
                is_diff(i,treedep)=true;
            end
        end
    end
end
end
