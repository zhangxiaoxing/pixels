function [permhat,permci]=binoperm(pos,tot,opt)
arguments
    pos (1,1) double
    tot (1,1) double {mustBeNonzero}
    opt.rpt (1,1) double = 100
end

permhat=pos./tot;
pool=zeros(tot,1);
pool(1:pos)=1;
permci=bootci(opt.rpt, @mean, pool);
end
