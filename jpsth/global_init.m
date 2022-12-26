global gather_config
gather_config=struct();
gather_config.fc_win=10;
gather_config.adjust_white_matter=true;
gather_config.corr_type='Pearson';
gather_config.fnsuffix='_10ms_adj_pearson';
gather_config.FC_thresh=100;

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% TODO:should be placed inside a structure
if strcmp(gather_config.corr_type,'Pearson')
    corr_ln_log='PearsonLinearLog';
    corr_log_log='PearsonLogLog';
else
    corr_ln_log='Spearman';
    corr_log_log='Spearman';
end