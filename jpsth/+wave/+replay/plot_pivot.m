function fh=plot_pivot(pivot_dict)
arguments
    pivot_dict = []
end

if isempty(pivot_dict)
    load(fullfile('binary','delay_vs_iti_pivot_spike.mat'),'pivot_dict');
end

out=cell2struct({[];[];[]},{'delay','iti','out_task'});
for fn=["delay","iti","out_task"]
    sucell=struct2cell(pivot_dict.(fn));
    for cii=1:numel(sucell)
        if ~isempty(sucell{cii}.values)
            out.(fn)=[out.(fn);mean(sucell{cii}.values)];
        end
    end
end

mm=[mean(out.delay),mean(out.iti),mean(out.out_task)];
stdd=[std(out.delay),std(out.iti),std(out.out_task)];
sem=stdd./sqrt([numel(out.delay),numel(out.iti),numel(out.out_task)]);

fh=figure('Position',[100,100,400,300]);
hold on;
bar(1:3,mm,'FaceColor','none')
errorbar(1:3,mm,sem,'k.','CapSize',12)
ylabel('Alternative path per spike');
xlim([0.25,3.75]);
set(gca,'XTick',1:3,'XTickLabel',{'Delay','ITI','out_task'},'XTickLabelRotation',90)

appendfig('tag','alternative path per spike, plot_pilot.m')
savefig(fh,fullfile('binary','delay_iti_alt_path_per_spk.fig'));
end

