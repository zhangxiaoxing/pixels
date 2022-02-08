% TODO cell assembly as algorithm
% Visualization

%% Rate selectivity index
% Selective
% SU

% FCSP

% Ring loops



%% Rate decoding
% Units vs decoding accuracy
% Selective congruent
% SU
su_data=bz.rings.cellasm.get_SU_data();

prefs=cellfun(@(x) mean(x('pref')),su_data);
nonps=cellfun(@(x) mean(x('nonp')),su_data);
freqsel=prefs>0.5;
%decoding
ringdec=bz.rings.cellasm.decode_one(su_data(freqsel),'unit_count',[10,50,100,500],'pre_pca',true);
cellfun(@(x) mean(x.result),ringdec)
cellfun(@(x) mean(x.su_count),ringdec)


% FCSP

% Ring loops
rings_data=bz.rings.cellasm.get_rings_data();

prefs=cellfun(@(x) mean(x('pref')),rings_data);
nonps=cellfun(@(x) mean(x('nonp')),rings_data);
freqsel=prefs>0.5;
%decoding
ringdec=bz.rings.cellasm.decode_one(rings_data(freqsel),'unit_count',[10,50,100,500],'pre_pca',true);
cellfun(@(x) mean(x.result),ringdec)
cellfun(@(x) mean(x.su_count),ringdec)

%selectivity index
selidx=(prefs(freqsel)-nonps(freqsel))./(prefs(freqsel)+nonps(freqsel));
%total number of su in rings
meta=cellfun(@(x) x('meta'),rings_data,'UniformOutput',false);
numel(unique(cell2mat(cellfun(@(x) x{4}*100000+x{5}, (meta(freqsel)).','UniformOutput',false))))



% Random-ML





function plot_one(data,shuf)
end


%% get data set, preferred vs non-preferred
% SU|feature:cell * label:map * trial:vector


function out=get_fcsp_data()
end



