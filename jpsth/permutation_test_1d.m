function p=permutation_test_1d(statsA,statsB,rpt)
     pool=[reshape(statsA,[],1);reshape(statsB,[],1)];
     numA=numel(statsA);
     nump=numel(pool);
     permDiff=nan(rpt,1);
     actualDiff=mean(statsA)-mean(statsB);
     for ii=1:rpt
        npool=randsample(pool,nump);
        permDiff(ii)=mean(npool(1:numA))-mean(npool((numA+1):end));
     end
     z=-abs(actualDiff-mean(permDiff))/std(permDiff);
     p=2*normcdf(z);
end