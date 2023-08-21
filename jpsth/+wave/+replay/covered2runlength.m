function [run_length,offset,onset]=covered2runlength(covered)
edges = find(diff([0;covered;0]==1));
onset = edges(1:2:end-1);  % Start indices
offset = edges(2:2:end);
run_length =(offset-onset)./10;  % Consecutive ones counts
end

