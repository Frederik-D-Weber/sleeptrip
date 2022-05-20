function [res_swsp_summary, res_swsp_channel_stat, res_nonswsp_channel_stat, res_nonspsw_channel_stat, res_swsp_event, res_nonswsp_event, res_nonspsw_event, res_excluded_sp_event, res_excluded_sw_event] = st_swsp(cfg, res_spindles_event, res_slowwaves_event)

% ST_SWSP find the slow-wave spindles (SW-spindles) from the resutls of st_spindles and st_slowwaves (e.g. ST_SPINDLES, ST_SLOWWAVES)
%
% Use as
%   [res_swsp_summary, res_swsp_channel_stat, res_nonswsp_channel_stat, res_nonspsw_channel_stat, res_swsp_event, res_nonswsp_event, res_nonspsw_event, res_excluded_sp_event, res_excluded_sw_event] = st_swsp(cfg, res_spindles_event, res_slowwaves_event)
%
% Optional configuration parameters are:
%
% cfg.SpindleTimePointColumn   = name of column with spindle time points in res_st_spindles.table (default = 'seconds_trough_max');
% cfg.SlowwaveTimePointColumn  = name of column with slow wave beginning time point in res_st_slowwaves.table (default = 'seconds_begin');
% cfg.SlowwaveTimeWindowOffsetTime  = offset in seconds from cfg.SlowwaveTimePointColumn (default = 0);
% cfg.SlowwaveTimePointColumn2 = name of column with slow wave ending time point in res_st_slowwaves.table (default = 'seconds_end');
% cfg.SlowwaveTimeWindowOffsetTime2  = offset in seconds from cfg.SlowwaveTimePointColumn2 (default = 0);
% cfg.groupby                  = cellstr with all the column names to group by in res_st_spindles and res_st_slowwaves (default = {'resnum', 'channel'}, e.g. = {'resnum', 'channel', stage_alt'} for also separating by sleep stages);
%
% See also ST_SPINDLES, ST_SLOWWAVES, ST_RES_FILTER

% Copyright (C) 2019-, Frederik D. Weber
%
% This file is part of SleepTrip, see http://www.sleeptrip.org
% for the documentation and details.
%
%    SleepTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    SleepTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    SleepTrip is a branch of FieldTrip, see http://www.fieldtriptoolbox.org
%    and adds funtionality to analyse sleep and polysomnographic data.
%    SleepTrip is under the same license conditions as FieldTrip.
%
%    You should have received a copy of the GNU General Public License
%    along with SleepTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
st_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set defaults
%cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.SpindleTimePointColumn  = ft_getopt(cfg, 'SpindleTimePointColumn', 'seconds_trough_max');
cfg.SlowwaveTimePointColumn  = ft_getopt(cfg, 'SlowwaveTimePointColumn', 'seconds_begin');
cfg.SlowwaveTimeWindowOffsetTime  = ft_getopt(cfg, 'SlowwaveTimeWindowOffsetTime', 0);
cfg.SlowwaveTimePointColumn2  = ft_getopt(cfg, 'SlowwaveTimePointColumn2', 'seconds_end');
cfg.SlowwaveTimeWindowOffsetTime2  = ft_getopt(cfg, 'SlowwaveTimeWindowOffsetTime2', 0);
cfg.groupby = ft_getopt(cfg, 'groupby',  {'resnum', 'channel'});

fprintf([functionname ' function initialized\n']);

cfg_co = cfg;
% the spindle time point is its maximal trough
cfg_co.EventsTestTimePointColumn = cfg.SpindleTimePointColumn;

% find SO spindles within the begin and end of slow waves, this
% assures that one spindle can target only one slow wave
cfg_co.UseSecondColumnAndOnlyOffsetsForTimeWindowTarget = 'yes';
cfg_co.EventsTargetTimePointColumn = cfg.SlowwaveTimePointColumn;
cfg_co.EventTargetTimeWindowOffsetTime = cfg.SlowwaveTimeWindowOffsetTime;
cfg_co.EventsTargetTimePointColumn2 = cfg.SlowwaveTimePointColumn2;
cfg_co.EventTargetTimeWindowOffsetTime2 = cfg.SlowwaveTimeWindowOffsetTime2;
cfg_co.column_prefix_test   =  'sp_';
cfg_co.column_prefix_target =  'sw_';

[res_swsp_summary, res_swsp_event, res_nonswsp_event, res_nonspsw_event, res_excluded_sp_event, res_excluded_sw_event] = st_coocc(cfg_co, res_spindles_event, res_slowwaves_event);

res_swsp_summary.ori = functionname;
res_swsp_event.ori = functionname;
res_nonswsp_event.ori = functionname;
res_nonspsw_event.ori = functionname;
res_excluded_sp_event.ori = functionname;
res_excluded_sw_event.ori = functionname;


res_swsp_summary.type = 'swsp_summary';
res_swsp_event.type = 'swsp_event';
res_nonswsp_event.type = 'nonswsp_event';
res_nonspsw_event.type = 'nonspsw_event';
res_excluded_sp_event.type = 'excluded_sp_event';
res_excluded_sw_event.type = 'excluded_sw_event';

res_swsp_event.table.delay_sp_time_minus_sw_offset = res_swsp_event.table.test_minus_target_delay;
res_swsp_event.table.test_minus_target_delay = [];

res_swsp_event.table.delay_sp_time_minus_sw_trough = res_swsp_event.table.([cfg_co.column_prefix_test cfg.SpindleTimePointColumn]) - res_swsp_event.table.sw_seconds_trough_max;

res_swsp_summary.table.sp_match_fraction = res_swsp_summary.table.sp_match_grouped./res_swsp_summary.table.sp_grouped;
res_swsp_summary.table.sw_match_fraction = res_swsp_summary.table.sp_match_grouped./res_swsp_summary.table.sw_ungrouped;


idcols = {};
for iGroupBy = 1:numel(cfg.groupby)
    idcols = cat(2,idcols,{[cfg_co.column_prefix_test cfg.groupby{iGroupBy}] [cfg_co.column_prefix_target cfg.groupby{iGroupBy}]});
end
%idcols = {'resnum','sp_channel','sw_channel'};
%idcols = {resnum','sp_channel','sw_channel','sp_stage_alt','sw_stage_alt'}

idcols_temp = idcols(ismember(idcols, res_swsp_event.table.Properties.VariableNames));

res_swsp_channel_stat = [];
res_swsp_channel_stat.ori = functionname;
res_swsp_channel_stat.type = 'swsp_channel_stat';
res_swsp_channel_stat.cfg = cfg;

datavars_temp = {...
    'sp_duration_seconds', 'sp_amplitude_peak2trough_max', 'sp_frequency_by_mean_pk_trgh_cnt_per_dur',...
    'sw_duration_seconds', 'sw_amplitude_peak2trough_max', 'sw_slope_to_trough_min_potential_per_second', 'sw_slope_zeroxing_to_trough_potential_per_second', 'sw_slope_trough_to_up_max_potential_per_second', 'sw_slope_trough_to_zeroxing_potential_per_second','sw_frequency_by_duration','sw_frequency_by_trough_to_peak_latency',...
    'delay_sp_time_minus_sw_offset','delay_sp_time_minus_sw_trough','sp_lengths_ROI_seconds'};

if size(res_swsp_event.table,1) > 0
    res_swsp_channel_stat.table = grpstats(res_swsp_event.table, idcols_temp,...
        {'mean','median','min','max','std'},...
        'DataVars',datavars_temp);
    res_swsp_channel_stat.table.Properties.VariableNames{(strcmp(res_swsp_channel_stat.table.Properties.VariableNames,'GroupCount'))} = 'count';
    res_swsp_channel_stat.table.density_per_minute = res_swsp_channel_stat.table.count ./ (res_swsp_channel_stat.table.mean_sp_lengths_ROI_seconds/60);
    res_swsp_channel_stat.table.Properties.RowNames = {};
else
    
    temptable = res_swsp_event.table(:,cat(2,idcols_temp,datavars_temp));
    temptab2 = cat(2,cell2table(idcols_temp),array2table(nan(1,size(datavars_temp,2))));
    %temptable(1,:) = temptab2(1,:);
    temptab2.Properties.VariableNames = temptable.Properties.VariableNames;
    res_swsp_channel_stat.table = grpstats(temptab2, idcols_temp,...
        {'mean','median','min','max','std'},...
        'DataVars',datavars_temp);
    
    res_swsp_channel_stat.table.Properties.VariableNames{(strcmp(res_swsp_channel_stat.table.Properties.VariableNames,'GroupCount'))} = 'count';
    res_swsp_channel_stat.table.density_per_minute = res_swsp_channel_stat.table.count ./ (res_swsp_channel_stat.table.mean_sp_lengths_ROI_seconds/60);
    res_swsp_channel_stat.table.Properties.RowNames = {};
    
    res_swsp_channel_stat.table(:,:) = [];
end


idcols_temp = idcols(ismember(idcols, res_nonswsp_event.table.Properties.VariableNames));
res_nonswsp_channel_stat = [];
res_nonswsp_channel_stat.ori = functionname;
res_nonswsp_channel_stat.type = 'nonswsp_channel_stat';
res_nonswsp_channel_stat.cft = cfg;

datavars_temp = {'sp_duration_seconds', 'sp_amplitude_peak2trough_max', 'sp_frequency_by_mean_pk_trgh_cnt_per_dur','sp_lengths_ROI_seconds'};

if size(res_nonswsp_event.table,1) > 0
    res_nonswsp_channel_stat.table = grpstats(res_nonswsp_event.table,idcols_temp,...
        {'mean','median','min','max','std'},...
        'DataVars',datavars_temp);
    res_nonswsp_channel_stat.table.Properties.VariableNames{(strcmp(res_nonswsp_channel_stat.table.Properties.VariableNames,'GroupCount'))} = 'count';
    res_nonswsp_channel_stat.table.density_per_minute = res_nonswsp_channel_stat.table.count ./ (res_nonswsp_channel_stat.table.mean_sp_lengths_ROI_seconds/60);
    res_nonswsp_channel_stat.table.Properties.RowNames = {};
else
    
    temptable = res_nonswsp_event.table(:,cat(2,idcols_temp,datavars_temp));
    temptab2 = cat(2,cell2table(idcols_temp),array2table(nan(1,size(datavars_temp,2))));
    %temptable(1,:) = temptab2(1,:);
    temptab2.Properties.VariableNames = temptable.Properties.VariableNames;
    
    res_nonswsp_channel_stat.table = grpstats(temptab2,idcols_temp,...
        {'mean','median','min','max','std'},...
        'DataVars',datavars_temp);
    res_nonswsp_channel_stat.table.Properties.VariableNames{(strcmp(res_nonswsp_channel_stat.table.Properties.VariableNames,'GroupCount'))} = 'count';
    res_nonswsp_channel_stat.table.density_per_minute = res_nonswsp_channel_stat.table.count ./ (res_nonswsp_channel_stat.table.mean_sp_lengths_ROI_seconds/60);
    res_nonswsp_channel_stat.table.Properties.RowNames = {};
    
    res_nonswsp_channel_stat.table(:,:) = [];
    
end

idcols_temp = idcols(ismember(idcols, res_nonspsw_event.table.Properties.VariableNames));
res_nonspsw_channel_stat = [];
res_nonspsw_channel_stat.ori = functionname;
res_nonspsw_channel_stat.type = 'nonspsw_channel_stat';
res_nonspsw_channel_stat.cft = cfg;

datavars_temp = {'sw_duration_seconds', 'sw_amplitude_peak2trough_max', 'sw_slope_to_trough_min_potential_per_second','sw_slope_zeroxing_to_trough_potential_per_second','sw_slope_trough_to_up_max_potential_per_second','sw_slope_trough_to_zeroxing_potential_per_second','sw_frequency_by_duration','sw_frequency_by_trough_to_peak_latency','sw_lengths_ROI_seconds'};

if size(res_nonspsw_event.table,1) > 0
    res_nonspsw_channel_stat.table = grpstats(res_nonspsw_event.table, idcols_temp,...
        {'mean','median','min','max','std'},...
        'DataVars',datavars_temp);
    res_nonspsw_channel_stat.table.Properties.VariableNames{(strcmp(res_nonspsw_channel_stat.table.Properties.VariableNames,'GroupCount'))} = 'count';
    res_nonspsw_channel_stat.table.density_per_minute = res_nonspsw_channel_stat.table.count ./ (res_nonspsw_channel_stat.table.mean_sw_lengths_ROI_seconds/60);
    res_nonspsw_channel_stat.table.Properties.RowNames = {};
else
    
    temptable = res_nonspsw_event.table(:,cat(2,idcols_temp,datavars_temp));
    temptab2 = cat(2,cell2table(idcols_temp),array2table(nan(1,size(datavars_temp,2))));
    %temptable(1,:) = temptab2(1,:);
    temptab2.Properties.VariableNames = temptable.Properties.VariableNames;
    
    res_nonspsw_channel_stat.table = grpstats(temptab2, idcols_temp,...
        {'mean','median','min','max','std'},...
        'DataVars',datavars_temp);
    
    res_nonspsw_channel_stat.table.Properties.VariableNames{(strcmp(res_nonspsw_channel_stat.table.Properties.VariableNames,'GroupCount'))} = 'count';
    res_nonspsw_channel_stat.table.density_per_minute = res_nonspsw_channel_stat.table.count ./ (res_nonspsw_channel_stat.table.mean_sw_lengths_ROI_seconds/60);
    res_nonspsw_channel_stat.table.Properties.RowNames = {};
    
    res_nonspsw_channel_stat.table(:,:) = [];
    
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
% ft_postamble previous data
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data



fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
