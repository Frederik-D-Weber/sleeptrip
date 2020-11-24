function [scoring res_events_normed] = st_scoringnormcycle(cfg, scoring, res_cycle, varargin)

% ST_SCORINGNORMCYCLE aligns a scoring/res_event structure like from
% ST_READ_SCORING/ST_SPINDLES/ST_SLOWAVES...
% with detected sleep cycles of this scoring as an output of
% ST_SLEEPCYCLES in a res_cycle structure to fit a new cycle scheme
% Use as
%   [scoring] = st_normcycle(cfg, scoring, res_cycle)
%   [scoring, res_events_normed] = st_scoringnormcycle(cfg, scoring, res_cycle, res_event)
%   [scoring, res_events_normed] = st_scoringnormcycle(cfg, scoring, res_cycle, res_event1, res_event2, ...)
%   [res_events_normed] = st_scoringnormcycle(cfg, scoring, res_cycle, res_event)
%   [res_events_normed] = st_scoringnormcycle(cfg, scoring, res_cycle, res_event1, res_event2, ...)
%
% Required configuration parameters are:
%   cfg.newcycledurations = either a string or a Cx1 vector giving the new
%                           cycle durations or a Cx2 matrix giving the NR and
%                           R cycle durations of each C cycle in the rows
%                           default is a 100x1 vector aligning to 90-min
%                           cycles
%
% if one or multiple res_event structures are provided the other optional
% arguments can be given:
%   cfg.eventdatatimecolumn = either a string or a cellstr of length in the
%                           number of res_event provided defining the 
%                           columnname with the event time in seconds 
%                           e.g. 'seconds_trough_max'
%                           default is any of the 'seconds_trough_max','te_seconds_trough_max', 'ta_seconds_trough_max', 'time', 'timepoint';
%
%   cfg.groupbyColumns =    either a string or a cell of cellstr of length in the
%                           number of cells given by res_event provided.
%                           Defines the columnnames for each res_event that
%                           the events shall be grouped by (e.g. the propertis of
%                           events)
%                           default = is {'channel'} for every res_event structure
%                           e.g. {'channel'} for a res_spindles_event structure
%                           or e.g. {{'channel','location'},{'channel'}}
%                           for an res_spindles_event and another res_event
%                           sturcture
%
%   cfg.meanColumns =       either a string or a cell of cellstr of length in the
%                           number of cells given by res_event provided.
%                           Defines the columnnames for each res_event that
%                           shall be averaged (e.g. the properties of
%                           events)
%                           e.g. {'duration_seconds','amplitude_peak2trough_max','frequency_by_mean_pk_trgh_cnt_per_dur'} for a res_spindles_event structure
%                           or e.g.
%                           {{'duration_seconds','amplitude_peak2trough_max'},{'duration_seconds','amplitude_peak2trough_max'}} for two res_spindles_event structures
%                           default = is {} for every res_event structure
%                           (i.e. no further properties but the count
%                           itself)
%
%  cfg.defaultvalues      = either a vector or a cell of vectors with the default values of the
%                           repective cfg.meanColumns, default =
%                           {nan(1,numel(cfg.meanColumns)), ...}
%  cfg.eventexcludeRNR    = cell of strings defining the exclusion of events that fall in the cycle parts of
%                           either R or NR for each event structure
%                           this happens BEFORE interpolation and smoothing
%                           example for 3 res_event sturctures
%                           cfg.excludeRNR = {'no', 'NR', 'R'}
%                           this will exclude the no events in the first
%                           res_event structure, all events in the NR cycle
%                           part of the second res_event structure, 
%                           and all events in the R cycle
%                           part of the third res_event structure
%                           default = {'no', ...}
%
%  cfg.eventmaskRNR       = cell of strings defining the masking of events that fall in the cycle parts of
%                           either R or NR for each event structure
%                           this happens AFTER interpolation and smoothing
%                           example for 3 res_event sturctures
%                           cfg.maskRNR = {'no', 'NR', 'R'}
%                           this will mask the no events in the first
%                           res_event structure, all events in the NR cycle
%                           part of the second res_event structure, 
%                           and all events in the R cycle
%                           part of the third res_event structure
%                           default = {'no', ...}
%
%  cfg.eventcycleinterpolmethd = interpolation method for res_event
%                           sturcture value-time normalization in cycle adjustments
%                           e.g. 'nearest' or 'linear' or 'pchip' see
%                           interp1 help for details on further methods
%                           default = 'nearest' 
%
%  cfg.eventepochinterpol = 'yes or 'no' if interpolation for res_event
%                           sturcture should be applied prior ot cycle
%                           interpolation to query a value for each epoch
%                           in the covered epoch range 
%                           default = 'no'
%  cfg.eventepochinterpolmethd = the method used for interpolation for res_event
%                           sturcture should be applied prior ot cycle
%                           interpolation to query a value for each epoch
%                           in the covered epoch range 
%                           e.g. 'nearest' or 'linear' or 'pchip' see
%                           interp1 help for details on further methods
%                           default = 'linear'
%
%  cfg.epochsmoother      = number of epochs being used to smooth values 
%                           from res_event prior to cycle interpolation;
%                           default = 1 (i.e. no smoothing)
%  cfg.epochsmoothermethd = method of epoch smoother
%                           e.g. 'moving' or 'lowess' or 'loess' 'rlowess' 'sgolay' see
%                           smooth function help for details on further methods
%                           default = 'moving';
%
%  cfg.resorientation =     string, orientation of res_events_normed epochs 
%                           either 'long' or 'wide' format default = 'wide'
%
%  cfg.considerexclusion =  string, either 'yes' or 'no' if the scoring.exclusion 
%                           should also be applied to events. default = 'yes';
%
%
% See also ST_SLEEPCYCLES, ST_READ_SCORING, ST_SCORINGCONVERT, ST_SLOWAVES,
% ST_SPINDLES

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

% FIXME check if the input res structures are from SleepTrip or result
% structures at all
% % check if the input data is valid for this function
% for i=1:length(varargin)
%   varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'raw', 'feedback', 'no', 'hassampleinfo', 'yes');
% end

% determine the dimensions of the data
Nres = length(varargin);

adjust_res = false;
if Nres>0
    adjust_res = true;
end

process_scoring = true;
if nargout == 1 && adjust_res
    process_scoring = false;
end



% set defaults
cfg.newcycledurations  = ft_getopt(cfg, 'newcycledurations', repmat(round(90*60/scoring.epochlength),100,1));

eventdatatimecolumn_candidates = {'seconds_trough_max', 'time', 'timepoint' ,'te_seconds_trough_max','test_seconds_trough_max','te_time','te_timepoint', 'ta_seconds_trough_max', 'target_seconds_trough_max', 'ta_time','ta_timepoint'};
if Nres>0
    if isfield(cfg,'eventdatatimecolumn')
        if (Nres ~= numel(cfg.eventdatatimecolumn))
            ft_error(['Number of column names in cfg.eventdatatimecolumn does not match the number of res_event structures = ' num2str(Nres) '\nCheck that there is a event column for each res_event sturcture given as an argument.'])
        end
    else
        ft_warning('attempting to automatically find cfg.eventdatatimecolumn:\n');
        cfg.eventdatatimecolumn = {};
        found = 0;
        for iResEvent = 1:Nres
            idxTimeColumns = ismember(eventdatatimecolumn_candidates,lower(varargin{iResEvent}.table.Properties.VariableNames));
            if any(idxTimeColumns)
                cfg.eventdatatimecolumn{iResEvent} = eventdatatimecolumn_candidates{find(idxTimeColumns,1)};
                ft_warning('found in res{%d} column ''%s''\n',iResEvent,cfg.eventdatatimecolumn{iResEvent});
                found = found + 1;
            end
        end
        if found ~= Nres
            ft_error('Could not find any matching columns that can be used for cfg.eventdatatimecolumn')
        end
    end
end

defaultGroupByColumns = true;
if Nres>0 
    if isfield(cfg,'groupbyColumns')
        if Nres > 1 && (Nres ~= numel(cfg.groupbyColumns))
            ft_error(['Number of cellstr in cfg.groupbyColumns does not match the number of res_event structures = ' num2str(Nres) '\nCheck that there is a cellstr for each res_event sturcture given as an argument.'])
        end
        defaultGroupByColumns = false;
    else
        
    cfg.groupbyColumns = {};
     for iResEvent = 1:Nres
         cfg.groupbyColumns{iResEvent} = {'channel'};
     end
    end
end

defaultMeanColumns = true;

if Nres>0 
    if isfield(cfg,'meanColumns')
        if Nres > 1 && (Nres ~= numel(cfg.meanColumns))
            ft_error(['Number of cellstr in cfg.meanColumns does not match the number of res_event structures = ' num2str(Nres) '\nCheck that there is a cellstr for each res_event sturcture given as an argument.'])
        end
        defaultMeanColumns = false;
    else
    cfg.meanColumns = {};
     for iResEvent = 1:Nres
         cfg.meanColumns{iResEvent} = {};
     end
    end
end


if Nres>0 
    if isfield(cfg,'defaultvalues')
        if Nres > 1 && (Nres ~= numel(cfg.defaultvalues))
            ft_error(['Number of cells in cfg.defaultvalues does not match the number of res_event structures = ' num2str(Nres) '\nCheck that there is a cell of values for each res_event sturcture given as an argument.'])
        end
    else
    cfg.defaultvalues = {};
     for iResEvent = 1:Nres
         cfg.defaultvalues{iResEvent} = nan(1,numel(cfg.meanColumns{iResEvent}));
     end
    end
end

if Nres>0 
    if isfield(cfg,'eventmaskRNR')
        if Nres > 1 && (Nres ~= numel(cfg.eventmaskRNR))
            ft_error(['Number of cells in cfg.eventmaskRNR does not match the number of res_event structures = ' num2str(Nres) '\nCheck that there is a cell of values for each res_event sturcture given as an argument.'])
        end
    else
    cfg.eventmaskRNR = {};
     for iResEvent = 1:Nres
         cfg.eventmaskRNR{iResEvent} = 'no';
     end
    end
end

if Nres>0 
    if isfield(cfg,'eventexcludeRNR')
        if Nres > 1 && (Nres ~= numel(cfg.eventexcludeRNR))
            ft_error(['Number of cells in cfg.eventexcludeRNR does not match the number of res_event structures = ' num2str(Nres) '\nCheck that there is a cell of values for each res_event sturcture given as an argument.'])
        end
    else
    cfg.eventexcludeRNR = {};
     for iResEvent = 1:Nres
         cfg.eventexcludeRNR{iResEvent} = 'no';
     end
    end
end

if adjust_res
    if ~iscell(cfg.eventdatatimecolumn)
        cfg.eventdatatimecolumn = {cfg.eventdatatimecolumn};
    end
    
    if ~iscell(cfg.groupbyColumns)
        cfg.groupbyColumns = {cfg.groupbyColumns};
    end
    
    if ~iscell(cfg.meanColumns)
        cfg.meanColumns = {cfg.meanColumns};
    end
    
    if ~iscell(cfg.defaultvalues)
        cfg.defaultvalues = {cfg.defaultvalues};
    end
    
    if (Nres == 1) && ~defaultGroupByColumns
        cfg.groupbyColumns = {cfg.groupbyColumns};
    end
    
    if (Nres == 1) && ~defaultMeanColumns
        cfg.meanColumns = {cfg.meanColumns};
    end
    

end

withinCycleAlign = false;
if size(cfg.newcycledurations,2) == 2
    withinCycleAlign = true;
end


cfg.eventcycleinterpolmethd = ft_getopt(cfg, 'eventcycleinterpolmethd', 'nearest');

cfg.eventepochinterpol = ft_getopt(cfg, 'eventepochinterpol', 'no');
cfg.eventepochinterpolmethd = ft_getopt(cfg, 'eventepochinterpolmethd', 'linear');
cfg.epochsmoother = ft_getopt(cfg, 'epochsmoother', 1);
cfg.epochsmoothermethd = ft_getopt(cfg, 'epochsmoothermethd', 'moving');


cfg.resorientation     = ft_getopt(cfg, 'resorientation', 'wide');

cfg.considerexclusion     = ft_getopt(cfg, 'considerexclusion', 'yes');



fprintf([functionname ' function initialized\n']);


completeCycleCount = sum(~isnan(res_cycle.table.endepoch));

if completeCycleCount > size(cfg.newcycledurations,1)
    ft_warning(['There are ' num2str(numel(res_cycle.table.cycle)) ' cycles defined in the res_cycle structure \nbut only the first ' num2str(size(cfg.newcycledurations,1)) ' cfg.newcycledurations will be used.|n Please check if the cycles match up.'])
    completeCycleCount = size(cfg.newcycledurations,1);
end


numberofepochsnew = nansum(nansum(cfg.newcycledurations(1:completeCycleCount,:),2));

if process_scoring
    
cfg_scc = [];
cfg_scc.to = 'number';


%scoring.label = unique(scoring.epochs)';
nLabels = numel(scoring.label);
scoring.prob = zeros(nLabels,numel(scoring.epochs));
for iLabel = 1:nLabels
    scoring.prob(iLabel,:) = strcmp(scoring.epochs,scoring.label{iLabel});
end

scoring_numbers = st_scoringconvert(cfg_scc,scoring);
scoring.numbers = cellfun(@str2num,scoring_numbers.epochs,'UniformOutput',1);


scorings_cycle_adjusted_number = [];
scorings_cycle_adjusted_excluded = [];
scorings_cycle_prop_adjusted = [];

event_values_cycle_adjusted_by_Result_and_Group_and_Column = {};
for iCycle = 1:completeCycleCount
    if withinCycleAlign
        iNR = (res_cycle.table.NRstartepoch(iCycle):res_cycle.table.NRendepoch(iCycle));
        iR = (res_cycle.table.Rstartepoch(iCycle):res_cycle.table.Rendepoch(iCycle));
        if isnan(iNR) % in case there is no NR and sleep started with a REM stage on sleep onset
            hyp_part_NR_adjusted = interp1_or_repeat(1,-1,linspace(min(iNR),max(iNR),cfg.newcycledurations(iCycle,1)),'nearest',NaN);
            hyp_part_NR_excluded_adjusted = interp1_or_repeat(1,single(logical(1)),linspace(min(iNR),max(iNR),cfg.newcycledurations(iCycle,1)),'nearest',NaN);
        else
            hyp_part_NR = scoring.numbers(iNR);
            hyp_part_NR_excluded = scoring.excluded(iNR);
            hyp_part_NR_adjusted = interp1_or_repeat(iNR,hyp_part_NR,linspace(min(iNR),max(iNR),cfg.newcycledurations(iCycle,1)),'nearest',NaN);
            hyp_part_NR_excluded_adjusted = interp1_or_repeat(iNR,single(hyp_part_NR_excluded),linspace(min(iNR),max(iNR),cfg.newcycledurations(iCycle,1)),'nearest',NaN);
        end
        
        if isnan(iR) % in case there is no NR and sleep started with a REM stage on sleep onset
            hyp_part_R_adjusted = interp1_or_repeat(1,-1,linspace(min(iR),max(iR),cfg.newcycledurations(iCycle,2)),'nearest',NaN);
            hyp_part_R_excluded_adjusted = interp1_or_repeat(1,single(logical(1)),linspace(min(iR),max(iR),cfg.newcycledurations(iCycle,2)),'nearest',NaN);
        else
            hyp_part_R = scoring.numbers(iR);
            hyp_part_R_excluded = scoring.excluded(iR);
            hyp_part_R_adjusted = interp1_or_repeat(iR,hyp_part_R,linspace(min(iR),max(iR),cfg.newcycledurations(iCycle,2)),'nearest',NaN);
            hyp_part_R_excluded_adjusted = interp1_or_repeat(iR,single(hyp_part_R_excluded),linspace(min(iR),max(iR),cfg.newcycledurations(iCycle,2)),'nearest',NaN);
        end
        
     
        
        chunk_scorings_NR_prop_adjusted = [];
        chunk_scorings_R_prop_adjusted = [];
        for iLabel = 1:nLabels
            if isnan(iNR) % in case there is no NR and sleep started with a REM stage on sleep onset
                hyp_part_NR_prob_adjusted = interp1_or_repeat(1,NaN,linspace(min(iNR),max(iNR),cfg.newcycledurations(iCycle,1)),'linear',NaN);
            else
                hyp_part_NR_prob = scoring.prob(iLabel,iNR);
                hyp_part_NR_prob_adjusted = interp1_or_repeat(iNR,hyp_part_NR_prob,linspace(min(iNR),max(iNR),cfg.newcycledurations(iCycle,1)),'linear',NaN);
            end
            
            if isnan(iR) % in case there is no NR and sleep started with a REM stage on sleep onset
                hyp_part_R_prob_adjusted  = interp1_or_repeat(1, NaN, linspace(min(iR), max(iR), cfg.newcycledurations(iCycle,2)),'linear',NaN);
            else
                hyp_part_R_prob = scoring.prob(iLabel,iR);
                hyp_part_R_prob_adjusted  = interp1_or_repeat(iR, hyp_part_R_prob, linspace(min(iR), max(iR), cfg.newcycledurations(iCycle,2)),'linear',NaN);
            end
    
            chunk_scorings_NR_prop_adjusted(iLabel,:) = hyp_part_NR_prob_adjusted;
            chunk_scorings_R_prop_adjusted(iLabel,:)  = hyp_part_R_prob_adjusted;

        end
        
        scorings_cycle_prop_adjusted   = cat(2,scorings_cycle_prop_adjusted,chunk_scorings_NR_prop_adjusted,chunk_scorings_R_prop_adjusted);
        scorings_cycle_adjusted_number = cat(2,scorings_cycle_adjusted_number,hyp_part_NR_adjusted,hyp_part_R_adjusted);
        scorings_cycle_adjusted_excluded = cat(2,scorings_cycle_adjusted_excluded,hyp_part_NR_excluded_adjusted,hyp_part_R_excluded_adjusted);
        
    else
        iC = (res_cycle.table.startepoch(iCycle):res_cycle.table.endepoch(iCycle));
        hyp_part_cycle = scoring.numbers(iC);
        hyp_part_cycle_excluded = scoring.excluded(iC);

        
        hyp_part_cycle_adjusted = interp1_or_repeat(iC,hyp_part_cycle,linspace(min(iC),max(iC),cfg.newcycledurations(iCycle,1)),'nearest',NaN);
        hyp_part_cycle_excluded_adjusted = interp1_or_repeat(iC,single(hyp_part_cycle_excluded),linspace(min(iC),max(iC),cfg.newcycledurations(iCycle,1)),'nearest',NaN);

        
        chunk_scorings_cycle_prop_adjusted = [];
        for iLabel = 1:nLabels
            hyp_part_cycle_prob = scoring.prob(iLabel,iC);
            hyp_part_cycle_prob_adjusted = interp1_or_repeat(iC,hyp_part_cycle_prob,linspace(min(iC),max(iC),cfg.newcycledurations(iCycle,1)),'linear',NaN);
            chunk_scorings_cycle_prop_adjusted(iLabel,:) = hyp_part_cycle_prob_adjusted;
        end
        
        scorings_cycle_prop_adjusted   = cat(2,scorings_cycle_prop_adjusted,chunk_scorings_cycle_prop_adjusted);
        scorings_cycle_adjusted_number = cat(2,scorings_cycle_adjusted_number,hyp_part_cycle_adjusted);
        scorings_cycle_adjusted_excluded = cat(2,scorings_cycle_adjusted_excluded,hyp_part_cycle_excluded_adjusted);

    end
end



scoringnew = scoring;
scoringnew.numbers = scorings_cycle_adjusted_number;
scoringnew.epochs = cellfun(@num2str,num2cell(scoringnew.numbers),'UniformOutput',0);
scoringnew.prob = scorings_cycle_prop_adjusted;
scoringnew.excluded = logical(scorings_cycle_adjusted_excluded);
scoringnew.standard = 'number';

end

res_events_normed = [];
if adjust_res
    epochNames = arrayfun(@(s) ['epoch_', num2str(s)], 1:numberofepochsnew, 'UniformOutput', false);
    tempIDnames = cat(2,{'reseventnum'},{'res_ori'},{'res_type'},{'groupby'},{'property'});
    tempvarnames = [tempIDnames,epochNames];
    res_scorings_cycle_normed_epochs_event_table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);

    %event_values_cycle_adjusted_by_Result_and_Group_and_Column = {};
    for iResEvent = 1:Nres
        maskRNR = cfg.eventmaskRNR{iResEvent};
        excludeRNR = cfg.eventexcludeRNR{iResEvent};
        
        res_event = varargin{iResEvent};
        
        res_event.table.epochnumber = floor((res_event.table.(cfg.eventdatatimecolumn{iResEvent})+scoring.dataoffset)/scoring.epochlength) + 1;
        
        if any(strcmp(res_event.table.Properties.VariableNames,'resnum'))
            groupBy = cat(2,{'resnum'},cfg.groupbyColumns{iResEvent},{'epochnumber'});
        else
            groupBy = cat(2,cfg.groupbyColumns{iResEvent},{'epochnumber'});
        end
        
        meanColumns = cfg.meanColumns{iResEvent};
        
        event_count = grpstats(res_event.table(:,cat(2,groupBy,meanColumns)),groupBy,'mean');
        event_count.Properties.VariableNames(strcmp(event_count.Properties.VariableNames,['GroupCount'])) = {'count'};
        event_count.Properties.VariableNames((end-numel(meanColumns)+1):end) = meanColumns;
        
        if istrue(cfg.considerexclusion)
            event_count(ismember(event_count.epochnumber,find(scoring.excluded)'),:) = [];
        end
        
             switch excludeRNR
                    case 'R'
                        idxRs = [];
                        for iCycle = 1:completeCycleCount
                            iR = (res_cycle.table.Rstartepoch(iCycle):res_cycle.table.Rendepoch(iCycle));
                            idxRs = [idxRs iR];
                        end
                        if ~isnan(idxRs)
                            event_count(ismember(event_count.epochnumber,idxRs'),:) = [];
                        end
                    case 'NR'
                        idxNRs = [];
                        for iCycle = 1:completeCycleCount
                            iNR = (res_cycle.table.NRstartepoch(iCycle):res_cycle.table.NRendepoch(iCycle));
                            idxNRs = [idxNRs iNR];
                        end
                        if ~isnan(idxNRs)
                            event_count(ismember(event_count.epochnumber,idxNRs'),:) = [];
                        end
                end
                
        completeCycleCount = sum(~isnan(res_cycle.table.endepoch));
        
        event_count_group = unique(event_count(:,cfg.groupbyColumns{iResEvent}));
        
        
        adjust_columns = cat(2,{'count'}, meanColumns);
        default_values = cat(2,0,cfg.defaultvalues{iResEvent});
        
        for iGroup = 1:size(event_count_group,1)
            
            groupValues = event_count_group(iGroup,:);
            
            matchIndicator = logical(ones(size(event_count,1),1));
            
            groupNames = event_count_group.Properties.VariableNames;
            for iComb = 1:numel(groupNames)
                %iComp = 1
                tempComp = groupNames{iComb};
                
                if iscell(event_count.(tempComp))
                    matchIndicator = matchIndicator & ( strcmp(event_count.(tempComp), groupValues.(tempComp)) );
                else
                    matchIndicator = matchIndicator & ( event_count.(tempComp) == groupValues.(tempComp) );
                end
            end
            
            event_count_group_column = event_count(matchIndicator,:);
            
            for iColumn = 1:numel(adjust_columns)
                adjust_column = adjust_columns{iColumn};
                default_value = default_values(iColumn);
                
                event_values_by_epoch = repmat(default_value,1,numel(scoring.epochs));
                if istrue(cfg.eventepochinterpol)
                	temp_epochnumbers = event_count_group_column.epochnumber;
                    temp_values = event_count_group_column.(adjust_column);
                    temp_epochnumbersquery = min(event_count_group_column.epochnumber):max(event_count_group_column.epochnumber);
                    temp_values_epochinterpol = interp1_or_repeat(temp_epochnumbers,temp_values,temp_epochnumbersquery,cfg.eventepochinterpolmethd,default_value);
                    event_values_by_epoch(temp_epochnumbersquery) = temp_values_epochinterpol;
                else
                    temp_epochnumbers = event_count_group_column.epochnumber;
                    temp_values = event_count_group_column.(adjust_column);
                    event_values_by_epoch(temp_epochnumbers) = temp_values;
                end
                
                if (cfg.epochsmoother > 1)
                    event_values_by_epoch = smooth(event_values_by_epoch,cfg.epochsmoother,cfg.epochsmoothermethd);
                end
                
                switch maskRNR
                    case 'R'
                        idxRs = [];
                        for iCycle = 1:completeCycleCount
                            iR = (res_cycle.table.Rstartepoch(iCycle):res_cycle.table.Rendepoch(iCycle));
                            idxRs = [idxRs iR];
                        end
                        if ~isnan(idxRs)
                            event_values_by_epoch(idxRs) = NaN;
                        end
                    case 'NR'
                        idxNRs = [];
                        for iCycle = 1:completeCycleCount
                            iNR = (res_cycle.table.NRstartepoch(iCycle):res_cycle.table.NRendepoch(iCycle));
                            idxNRs = [idxNRs iNR];
                        end
                        if ~isnan(idxNRs)
                            event_values_by_epoch(idxNRs) = NaN;
                        end
                end
                
                event_values_cycle_adjusted = [];
                for iCycle = 1:completeCycleCount
                    if withinCycleAlign
                        iNR = (res_cycle.table.NRstartepoch(iCycle):res_cycle.table.NRendepoch(iCycle));
                        iR = (res_cycle.table.Rstartepoch(iCycle):res_cycle.table.Rendepoch(iCycle));
                        
                        if isnan(iNR)
                            event_part_NR_adjusted = interp1_or_repeat(1,default_value,linspace(min(iNR),max(iNR),cfg.newcycledurations(iCycle,1)),cfg.eventcycleinterpolmethd,default_value);
                        else
                            event_part_NR = event_values_by_epoch(iNR);
                            event_part_NR_adjusted = interp1_or_repeat(iNR,event_part_NR,linspace(min(iNR),max(iNR),cfg.newcycledurations(iCycle,1)),cfg.eventcycleinterpolmethd,default_value);
                        end
                        
                        if isnan(iR)
                            event_part_R_adjusted = interp1_or_repeat(1,default_value,linspace(min(iR),max(iR),cfg.newcycledurations(iCycle,2)),cfg.eventcycleinterpolmethd,default_value);
                        else
                            event_part_R = event_values_by_epoch(iR);
                            event_part_R_adjusted = interp1_or_repeat(iR,event_part_R,linspace(min(iR),max(iR),cfg.newcycledurations(iCycle,2)),cfg.eventcycleinterpolmethd,default_value);
                        end
                        
                        event_values_cycle_adjusted = cat(2,event_values_cycle_adjusted,event_part_NR_adjusted,event_part_R_adjusted);
                    else
                        iC = (res_cycle.table.startepoch(iCycle):res_cycle.table.endepoch(iCycle));
                        event_part_cycle = event_values_by_epoch(iC);
                        
                        event_part_cycle_adjusted = interp1_or_repeat(iC,event_part_cycle,linspace(min(iC),max(iC),cfg.newcycledurations(iCycle,1)),cfg.eventcycleinterpolmethd,default_value);
                        event_values_cycle_adjusted = cat(2,event_values_cycle_adjusted,event_part_cycle_adjusted);
                    end
                end
                
                
        
                
  
        groupNameJoined = strjoin(cellfun(@num2str,table2cell(groupValues),'UniformOutput',false),' + ');
        res_scorings_cycle_normed_epochs_event_table = ...
        cat(1,res_scorings_cycle_normed_epochs_event_table,...
              cat(2,table(iResEvent,{res_event.ori},{res_event.type},{groupNameJoined},{adjust_column},'VariableNames',tempIDnames),...
                    array2table(event_values_cycle_adjusted,'VariableNames',epochNames)));

                
                
                
                
                
                
                %event_values_cycle_adjusted_by_Result_and_Group_and_Column{iResEvent,iGroup,iColumn} = event_values_cycle_adjusted;
            end
        end
        
        
    end
    
if strcmp(cfg.resorientation,'long')
res_scorings_cycle_normed_epochs_event_table = stack(res_scorings_cycle_normed_epochs_event_table,6:size(res_scorings_cycle_normed_epochs_event_table,2),...
          'NewDataVariableName','value',...
          'IndexVariableName','epoch');
end

res_events_normed = [];
res_events_normed.ori = functionname;
res_events_normed.type = ['res_events_normed_', cfg.resorientation];
res_events_normed.cfg = cfg;
res_events_normed.table = res_scorings_cycle_normed_epochs_event_table;

end


if process_scoring
    cfg = [];
    cfg.to = scoring.standard;
    scoring = st_scoringconvert(cfg, scoringnew);
else
    scoring = res_events_normed;
end




fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end

function vq = interp1_or_repeat(x,v,xq,method,defaultvalue) 

 if numel(x) > 1
    vq = interp1(x,v,xq,method);
 else
    if numel(x) == 1
        vq = repmat(v,size(xq));
    else
        vq = repmat(defaultvalue,size(xq));
    end
 end
end
