function T_simple=st_scoringdescriptives_select(cfg)

ft_checkconfig(cfg,'required',{'sleepArchitecture'});
sleepArchitecture=cfg.sleepArchitecture;

cfg.desiredSleepVariables  = ft_getopt(cfg, 'desiredSleepVariables', []);
desiredSleepVars=cfg.desiredSleepVariables;

metaVars={...
    'scoring_duration_min',...
    'sleep_opportunity_on_exists',...
    'sleep_opportunity_off_exists',...
    'sleep_opportunity_on_min',...
    'sleep_opportunity_off_min',...
    'sleep_onset_time_after_scoring_start_min',...
    'sleep_offset_time_after_scoring_start_min'};

basicVars={...
    'sleep_opportunity_on_to_off_min',...
    'total_sleep_period_duration_min',...
    'total_sleep_duration_in_sleep_period_min',...
    'sleep_eff_total_sleep_dur_in_sleep_prd_perc_of_sleep_opport',...
    'sleep_eff_total_sleep_dur_in_sleep_prd_perc_of_sleep_prd',...
    'sleep_onset_delay_min',...
    'N1_of_sleep_period_min',...
    'N2_of_sleep_period_min',...
    'N3_of_sleep_period_min',...
    'R_of_sleep_period_min',...
    'Wake_after_sleep_onset_of_sleep_period_min',...
    'N1_perc_of_sleep_period',...
    'N2_perc_of_sleep_period',...
    'N3_perc_of_sleep_period',...
    'R_perc_of_sleep_period',...
    'Wake_after_sleep_onset_perc_of_sleep_period',...
    'N1_delay_min',...
    'N2_delay_min',...
    'SWS_delay_min',...
    'R_delay_min'};

    arousalVars={...
    'arousals_density_per_min_sleep_in_sleep_period',...
    'arousals_density_per_min_N1_in_sleep_period',...
    'arousals_density_per_min_N2_in_sleep_period',...
    'arousals_density_per_min_N3_in_sleep_period',...
    'arousals_density_per_min_R_in_sleep_period'};


if isempty(desiredSleepVars)

    desiredSleepVars=[basicVars arousalVars];

end


T_ori=sleepArchitecture.table;
T_simple=T_ori(1,desiredSleepVars);