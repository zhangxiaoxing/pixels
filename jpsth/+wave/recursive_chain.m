function out=recursive_chain(chain,prior)
out=struct();
[out.cids,out.tcoms]=deal([]);
if numel(chain)>1
    if isempty(prior)
        for ii=1:size(chain{1}{1},1)
            stage.cids={chain{1}{1}(ii,:)};
            stage.tcoms={chain{1}{2}(ii,:)};
            rtn=wave.recursive_chain(chain(2:end),stage);
            out.cids=[out.cids;rtn.cids];
            out.tcoms=[out.tcoms;rtn.tcoms];
        end
    else
        cutoff=true;
        for ii=1:size(chain{1}{1},1)
            if chain{1}{1}(ii,1)==prior.cids{1}(end)
                cutoff=false;
                stage.cids={[prior.cids{1},chain{1}{1}(ii,2)]};
                stage.tcoms={[prior.tcoms{1},chain{1}{2}(ii,2)]};
                rtn=wave.recursive_chain(chain(2:end),stage);
                out.cids=[out.cids;rtn.cids];
                out.tcoms=[out.tcoms;rtn.tcoms];
            end
        end
        if cutoff
            if numel(prior.cids{1})>2
                out.cids=[out.cids;prior.cids];
                out.tcoms=[out.tcoms;prior.tcoms];
            end
        end
    end
else
    cutoff=true;
    for ii=1:size(chain{1}{1},1)
        if chain{1}{1}(ii,1)==prior.cids{1}(end)
            cutoff=false;
            out.cids={[prior.cids{1},chain{1}{1}(ii,2)]};
            out.tcoms={[prior.tcoms{1},chain{1}{2}(ii,2)]};
        end
    end
    if cutoff
        if numel(prior.cids{1})>2
            out.cids=[out.cids;prior.cids];
            out.tcoms=[out.tcoms;prior.tcoms];
        end
    end
end