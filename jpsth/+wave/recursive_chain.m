function out=recursive_chain(chain,prior)
out=[];
if numel(chain)>1
    if isempty(prior)
        for ii=1:size(chain{1}{1})
            out=[out;wave.recursive_chain(chain(2:end),[chain{1}{1}(ii,:);chain{1}{2}(ii,:)])];
        end
    else
        for ii=1:size(chain{1}{1})
            if chain{1}{1}(ii,1)==prior(1,end)
                out=[out;wave.recursive_chain(chain(2:end),[prior(1,:),chain{1}{1}(ii,2);prior(2,:),chain{1}{2}(ii,2)])];
            end
        end
    end
else
    for ii=1:size(chain{1}{1})
        if chain{1}{1}(ii,1)==prior(1,end)
            out=[out;[prior(1,:),chain{1}{1}(ii,2);prior(2,:),chain{1}{2}(ii,2)]];
        end
    end
end