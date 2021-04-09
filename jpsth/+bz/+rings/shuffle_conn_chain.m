function outconn=shuffle_conn_chain(in)
outconn=nan(size(in));
if nnz(sel)>0
    shufsel=randperm(size(in,1));
    shufdata=inpair(selpair(shufsel(1:nnz(sel))),:);
    flipsel=randi(2,size(shufdata,1),1)>1;
    shufdata(flipsel,:)=shufdata(flipsel,[2,1]);
    outconn(sel,:)=shufdata;
end
end