function plot_pivot(pivot_dict)
out=cell2struct({[];[];[]},{'delay','iti','surround'});
for fn=["delay","iti","surround"]
    sucell=struct2cell(pivot_dict.(fn));
    for cii=1:numel(sucell)
        if ~isempty(sucell{cii}.values)
            out.(fn)=[out.(fn);mean(sucell{cii}.values)];
        end
    end
end

mm=[mean(out.delay),mean(out.iti),mean(out.surround)];
stdd=[std(out.delay),std(out.iti),std(out.surround)];
sem=stdd./sqrt([numel(out.delay),numel(out.iti),numel(out.surround)]);

fh=figure('Position',[100,100,400,300])
hold on;
bar(1:3,mm,'FaceColor','none')
errorbar(1:3,mm,sem,'k.','CapSize',12)
ylabel('Alternative path per spike');
xlim([0.25,3.75]);
set(gca,'XTick',1:3,'XTickLabel',{'Delay','ITI','Surround'},'XTickLabelRotation',90)

end

