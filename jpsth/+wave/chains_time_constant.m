load('chains_mix.mat','chains_uf');
chains=chains_uf;
clear chains_uf;
wids=reshape(unique(chains.wave),1,[]);
t_span=struct();
for dur=3:3:6
    for wid=wids
        lsel=cellfun(@(x) numel(x)>4,chains.cids);
        wsel=chains.wave==wid & chains.dur==dur;
        dtcom=cellfun(@(x) x(end)-x(1),chains.tcoms(lsel & wsel));
%         if ~isfield(t_span,wid+num2str(dur))
%             t_span.(wid+num2str(dur))=[];
%         end
        if ~isempty(dtcom)
            t_span.("d"+num2str(dur)).(wid)=dtcom;
        end
    end
end

xx=[0:200:1000,2000:1000:6000];
hist3=histcounts([t_span.d3.olf_s1;t_span.d3.olf_s2;t_span.d3.s1d3;t_span.d3.s2d3]*1000/4,xx,'Normalization','pdf');
hist6=histcounts([t_span.d6.olf_s1;t_span.d6.olf_s2;t_span.d6.s1d6;t_span.d6.s2d6]*1000/4,xx,'Normalization','pdf');

figure()
hold on
plot([100:200:900,1500:1000:5500],hist3,'--k')
plot([100:200:900,1500:1000:5500],hist6,'-k')
set(gca(),'XScale','log','YScale','log');
ylim([1e-5,0.1])