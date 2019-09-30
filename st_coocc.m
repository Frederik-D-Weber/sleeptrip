function [res_summary, res_match_test_target, res_mismatch_test, res_mismatch_target, res_excluded_test, res_excluded_target] = st_coocc(cfg, res_test, res_target)

% ST_COOCC measures the co-occurrence of test within target events
% described by their time points, identifiers and properties.
% Test events are checked if they fall within a defined timewindow arround
% target events and thus either match them or mismatch them.
% It is based on output of event functions like ST_SPINDLES and
% ST_SLOWWAVES, and can be used to detect time-locking of e.g. spindles and
% slow waves, or prevalence of those events in multiple channels. It forms
% the basis for traveling wave/event analysis.
%
% The function takes two result files one test and one target one (the order or argument is
% very important).
%
% Use as
%   [res_summary, res_match_test_target, res_mismatch_test, res_mismatch_target, res_excluded_test, res_excluded_target] = st_coocc(cfg, res_test, res_target)
%   [res_summary, res_match_test_target, res_mismatch_test, res_mismatch_target] = st_coocc(cfg, res_test, res_target)
%   [res_summary, res_match_test_target, res_mismatch_test] = st_coocc(cfg, res_test, res_target)
%   [res_summary, res_match_test_target] = st_coocc(cfg, res_test, res_target)
%   [res_summary] = st_coocc(cfg, res_test, res_target)
%
% Required configuration parameters are:
%
%   cfg.EventsTestTimePointColumn = column name for event value in test
%                                 event result table of res_test
%                                 e.g. 'seconds_trough_max'
%   cfg.EventsTargetTimePointColumn = column name for event value in target
%                                 event result table of res_target
%                                 e.g. 'seconds_trough_max'
%                                 alternatively something like
%                                'seconds_begin' (for example cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow = 'yes')
%
% Optional configuration parameters are:
%   cfg.EventsTestCompareColumns  = either a string or a Nx1 cell-array with strings
%                                that indicate the columns used for matching test
%                                events against events columns in
%                                target events files,
%                                must match the size of cfg.EventsTargetCompareColumns
%                                (default = {'resnum', 'channel'}),
%  cfg.EventsTargetCompareColumns  = either a string or a Nx1 cell-array with strings
%                                that indicate the columns used for matching test
%                                events against events columns in
%                                target events files,
%                                must match the size of cfg.EventsTargetCompareColumns
%                               (default = {'resnum', 'channel'})
%  cfg.EventsTestGroupSummaryByColumns = column names for grouping the test event stats in the summary
%                               (default = {'resnum', 'channel', 'stage_alt'})
%  cfg.EventTargetTimeWindowOffsetTime = time points offset of the time window test events can fall into
%                               with respect to target events time points in units used for
%                               cfg.EventsTargetCompareColumns (default = 0)
%  cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow = if a second column and a second offset is used to
%                               define the time window around target events,
%                               either 'yes' or 'no' (default = 'no')
%                               if 'yes' the cfg.EventsTestTimePointColumn and cfg.EventsTestTimePointColumn2
%                                 and respective cfg.EventTargetTimeWindowOffsetTime and cfg.EventTargetTimeWindowOffsetTime2
%                                 are used to define the left and right boundaries for the time window
%                               if 'no' then instead cfg.EventsTestTimePointColumn and cfg.EventTargetTimeWindowOffsetTime
%                                 and cfg.EventTargetTimeWindowPreOffsetTime and cfg.EventTargetTimeWindowPostOffsetTime
% cfg.EventsTargetTimePointColumn2 = column name for event value in target
%                               event result table of res_target
%                               e.g. 'seconds_trough_max' or cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow is set to yes
% cfg.EventTargetTimeWindowOffsetTime2 = the second time point offset of the time window test events can fall into with respect to target events.
%                               The time points are in in units used for cfg.EventsTargetCompareColumns
%                               and is only used if cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow
%                               is set to 'yes' (default = 0)
% cfg.EventTargetTimeWindowPreOffsetTime = time (in seconds) before target event plus offset
%                               to define time window that test events can fall into
%                               with respect to target events.
%                               The time points are in units used for cfg.EventsTargetCompareColumns
%                               (default = 1)
% cfg.EventTargetTimeWindowPostOffsetTime = time (in seconds) after target event plus offset
%                               to define time window that test events can fall into
%                               with respect to target events.
%                               The time points in units used for cfg.EventsTargetCompareColumns
%                               (default = 1)
% cfg.EventTestTimePointOffsetTime = time point offset (in seconds) of test events
%                               in units used for cfg.EventsTestTimePointColumn
%                               (default = 0)
%
% cfg.MismatchIdenticalEvents = string indicating if events should be not matched
%                               when they are identical according to a comparison
%                               of the columns defined in
%                               cfg.EventsTestIDColumns vs cfg.EventsTargetIDColumns.
%                               This comes in handy for searching events between channels.
%                               either 'yes' or 'no' (default = 'no')
% cfg.MismatchDuplicateTestTargetOrTargetTestMatchingEvents = string indicating
%                               if duplicate matching test and target events according to
%                               cfg.EventsTestIDColumns vs cfg.EventsTargetIDColumns
%                               should be removed from results, note that the
%                               cfg.EventsTestIDColumns and cfg.EventsTargetIDColumns
%                               do NOT have to match cfg.EventsTestCompareColumns and cfg.EventsTargetCompareColumns.
%                               either 'yes' or 'no' (default = 'no')
% cfg.EventsTestIDColumns     = either a string or a Nx1 cell-array with strings
%                               that indicate the columns in res_test input used for identifying events, i.e. for mismatching identical or removing duplicate matches test
%                               (default = {'resnum', 'channel', cfg.EventsTestTimePointColumn})
% cfg.EventsTargetIDColumns   = either a string or a Nx1 cell-array with strings
%                               that indicate the columns in res_target input used for identifying events, i.e. for mismatching identical or removing duplicate matches test
%                               (default = {'resnum', 'channel', cfg.EventsTargetTimePointColumn})
%  cfg.EventsTestFilterForColumns = either a string or cellstr of the column(s) you want to filter
%                               in the test event (res_test) for, e.g. {'resnum', 'channel'}
%  cfg.EventsTestFilterValues = either a string/value or a cell of cells/cellstr,
%                               with the corresponding values to
%                               the columns defined cfg.EventsTestFilterForColumns,
%                               e.g. {{1}, {'C3', 'C4'}}
%  cfg.EventsTargetFilterForColumns = either a string/value or cellstr of the column(s) you want to filter
%                               in the target event (res_target) for, e.g. {'resnum', 'channel'}
%  cfg.EventsTargetFilterValues = either a string or a cell of cells/cellstr, with the corresponding values to
%                               the columns defined cfg.EventsTargetFilterForColumns,
%                               e.g. {{1}, {'C3', 'C4'}}
%  cfg.concatbuffersize       = size for concatenation buffer of results lets the datasets not get to
%                               big for concatenation of result this give the maximal
%                               number of lines a dataset can have for concatenation
%                               and when low can increase performance for >10000 events.
%                               (default = 100)
%
%
% See also ST_SPINDLES, ST_SLOWWAVES

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
memtic
st = dbstack;
functionname = st.name;
fprintf([functionname ' function started\n']);


% set defaults
cfg.EventsTestCompareColumns  = ft_getopt(cfg, 'EventsTestCompareColumns', {'resnum', 'channel'});
cfg.EventsTargetCompareColumns  = ft_getopt(cfg, 'EventsTargetCompareColumns', {'resnum', 'channel'});

cfg.EventsTestGroupSummaryByColumns  = ft_getopt(cfg, 'EventsTestGroupSummaryByColumns', {'resnum', 'channel', 'stage_alt'});
cfg.EventTargetTimeWindowOffsetTime = ft_getopt(cfg, 'EventTargetTimeWindowOffsetTime', 0);

cfg.EventsTestFilterForColumns  = ft_getopt(cfg, 'EventsTestFilterForColumns', {});
cfg.EventsTestFilterValues  = ft_getopt(cfg, 'EventsTestFilterValues', {});
cfg.EventsTargetFilterForColumns  = ft_getopt(cfg, 'EventsTargetFilterForColumns', {});
cfg.EventsTargetFilterValues  = ft_getopt(cfg, 'EventsTargetFilterValues', {});

cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow  = ft_getopt(cfg, 'UseSecondColumnAndOnlyOffsetsForTimeWindow', 'no');
cfg.EventsTargetTimePointColumn2  = ft_getopt(cfg, 'EventsTargetTimePointColumn2', 'seconds_end');
cfg.EventTargetTimeWindowOffsetTime2  = ft_getopt(cfg, 'EventTargetTimeWindowOffsetTime2', 0);
cfg.EventTargetTimeWindowPreOffsetTime  = ft_getopt(cfg, 'EventTargetTimeWindowPreOffsetTime', 1);
cfg.EventTargetTimeWindowPostOffsetTime  = ft_getopt(cfg, 'EventTargetTimeWindowPostOffsetTime', 1);
cfg.EventTestTimePointOffsetTime  = ft_getopt(cfg, 'EventTestTimePointOffsetTime', 0);

cfg.MismatchIdenticalEvents  = ft_getopt(cfg, 'MismatchIdenticalEvents', 'no');
cfg.MismatchDuplicateTestTargetOrTargetTestMatchingEvents  = ft_getopt(cfg, 'MismatchDuplicateTestTargetOrTargetTestMatchingEvents', 'no');
cfg.EventsTestIDColumns  = ft_getopt(cfg, 'EventsTestIDColumns', {'resnum', 'channel' cfg.EventsTestTimePointColumn});
if istrue(cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow)
    cfg.EventsTargetIDColumns  = ft_getopt(cfg, 'EventsTargetIDColumns', {'resnum', 'channel' cfg.EventsTargetTimePointColumn cfg.EventsTargetTimePointColumn2});
else
    cfg.EventsTargetIDColumns  = ft_getopt(cfg, 'EventsTargetIDColumns', {'resnum', 'channel' cfg.EventsTargetTimePointColumn});
end


cfg.concatbuffersize  = ft_getopt(cfg, 'concatbuffersize', 100);





if istrue(cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow) && istrue(cfg.MismatchIdenticalEvents)
    ft_error(['IMPORTANT: cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow = ''ýes'' and cfg.MismatchIdenticalEvents = ''ýes'' is not supported.'  ])
end


if ~iscell(cfg.EventsTestCompareColumns)
    cfg.EventsTestCompareColumns = {cfg.EventsTestCompareColumns};
end

if ~iscell(cfg.EventsTargetCompareColumns)
    cfg.EventsTargetCompareColumns = {cfg.EventsTargetCompareColumns};
end

if ~iscell(cfg.EventsTestFilterForColumns)
    cfg.EventsTestFilterForColumns = {cfg.EventsTestFilterForColumns};
end

if ~iscell(cfg.EventsTestFilterValues)
    cfg.EventsTestFilterValues = {cfg.EventsTestFilterValues};
end

if ~iscell(cfg.EventsTargetFilterValues)
    cfg.EventsTargetFilterValues = {cfg.EventsTargetFilterValues};
end

tempidx = ismember(cfg.EventsTestCompareColumns, res_test.table.Properties.VariableNames);
cfg.EventsTestCompareColumns = cfg.EventsTestCompareColumns(tempidx);
if ~all(tempidx)
    ft_warning(['dropped not existent column with the names ' strjoin(cfg.EventsTestCompareColumns{~tempidx},' ') 'in cfg.EventsTestCompareColumns']);
end

tempidx = ismember(cfg.EventsTargetCompareColumns, res_target.table.Properties.VariableNames);
cfg.EventsTestCompareColumns = cfg.EventsTestCompareColumns(tempidx);
if ~all(tempidx)
    ft_warning(['dropped not existent column with the names ' strjoin(cfg.EventsTargetCompareColumns{~tempidx},' ') 'in cfg.EventsTestCompareColumns']);
end

if ~(all(size(cfg.EventsTestCompareColumns) == size(cfg.EventsTargetCompareColumns)))
    error('number of test events and target events columns do not aggree')
end

if isfield(cfg,'EventsTestFilterForColumns')
    if ~isfield(cfg,'EventsTestFilterValues')
        ft_error(['if cfg.EventsTestFilterForColumns is defined then cfg.EventsTestFilterValues also needs to be defined'])
    end
    if numel(cfg.EventsTestFilterForColumns) ~= numel(cfg.EventsTestFilterValues)
        ft_error(['cfg.EventsTestFilterForColumns needs to have the same amount of elements as cfg.EventsTestFilterValues'])
    end
    
    %     for iCol = 1:numel(cfg.EventsTestFilterForColumns)
    %         fv = res_test.table{:,{cfg.EventsTestFilterForColumns{iCol}}};
    %         if ~isempty(fv)
    %             res_test.table = res_test.table{ismember(fv,cfg.EventsTestFilterValues{iCol}),:};
    %         end
    %     end
end

if isfield(cfg,'EventsTargetFilterForColumns')
    if ~isfield(cfg,'EventsTargetFilterValues')
        ft_error(['if cfg.EventsTargetFilterForColumns is defined then cfg.EventsTargetFilterValues also needs to be defined'])
    end
    if numel(cfg.EventsTargetFilterForColumns) ~= numel(cfg.EventsTargetFilterValues)
        ft_error(['cfg.EventsTargetFilterForColumns needs to have the same amount of elements as cfg.EventsTargetFilterValues'])
    end
    
    %     for iCol = 1:numel(cfg.EventsTargetFilterForColumns)
    %         fv = res_target.table{:,{cfg.EventsTargetFilterForColumns{iCol}}};
    %         if ~isempty(fv)
    %             res_target.table = res_target.table{ismember(fv,cfg.EventsTargetFilterValues{iCol}),:};
    %         end
    %     end
end

tempidx = ismember(cfg.EventsTestFilterForColumns, res_test.table.Properties.VariableNames);
cfg.EventsTestFilterForColumns = cfg.EventsTestFilterForColumns(tempidx);
cfg.EventsTestFilterValues = cfg.EventsTestFilterValues(tempidx);
if ~all(tempidx)
    ft_warning(['dropped not existent column with the names ' strjoin(cfg.EventsTestFilterForColumns{~tempidx},' ') 'in cfg.EventsTestFilterForColumns']);
end

tempidx = ismember(cfg.EventsTargetFilterForColumns, res_target.table.Properties.VariableNames);
cfg.EventsTargetFilterForColumns = cfg.EventsTargetFilterForColumns(tempidx);
cfg.EventsTargetFilterValues = cfg.EventsTargetFilterValues(tempidx);
if ~all(tempidx)
    ft_warning(['dropped not existent column with the names ' strjoin(cfg.EventsTargetFilterForColumns{~tempidx},' ') 'in cfg.EventsTargetFilterForColumns']);
end

tempidx = ismember(cfg.EventsTestIDColumns, res_test.table.Properties.VariableNames);
cfg.EventsTestIDColumns = cfg.EventsTestIDColumns(tempidx);
if ~all(tempidx)
    ft_warning(['dropped not existent column with the names ' strjoin(cfg.EventsTestIDColumns{~tempidx},' ') 'in cfg.EventsTestIDColumns']);
end

tempidx = ismember(cfg.EventsTargetIDColumns, res_target.table.Properties.VariableNames);
cfg.EventsTargetIDColumns = cfg.EventsTargetIDColumns(tempidx);
if ~all(tempidx)
    ft_warning(['dropped not existent column with the names ' strjoin(cfg.EventsTargetIDColumns{~tempidx},' ') 'in cfg.EventsTargetIDColumns']);
end

if (isfield(cfg,'EventsTestIDColumns') && ~isfield(cfg,'EventsTargetIDColumns')) || (~isfield(cfg,'EventsTestIDColumns') && isfield(cfg,'EventsTargetIDColumns'))
    error('cfg.EventsTestIDColumns and cfg.EventsTargetIDColumns need to be defined together')
end

if isfield(cfg,'EventsTestIDColumns') && isfield(cfg,'EventsTargetIDColumns') && ~istrue(cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow)
    if ~(all(size(cfg.EventsTestIDColumns) == size(cfg.EventsTargetIDColumns)))
        error('number of columns in cfg.EventsTestIDColumns and cfg.EventsTargetIDColumns do not aggree')
    end
end

fprintf([functionname ' function initialized\n']);

GroupByConcatString = '<#%%#>';


eventsTestIDname   = [res_test.ori '_' res_test.type];
eventsTargetIDname = [res_target.ori '_' res_target.type];

%TODO check if takes more memory than direct reference
EventsTestTable = res_test.table;
EventsTargetTable = res_target.table;

matchIndicator = ones(size(EventsTestTable,1),1);
if ~isempty(cfg.EventsTestFilterForColumns)
    matchIndicator = zeros(size(EventsTestTable,1),1);
    for iCol = 1:numel(cfg.EventsTestFilterForColumns)
        tempEventsTestFilterForColumn = cfg.EventsTestFilterForColumns{iCol};
        for iComb = 1:numel(cfg.EventsTestFilterValues)
            %iComp = 1
            tempCompTest = cfg.EventsTestFilterValues{iComb};
            %tempCompTarget = cfg.EventsTargetFilterValues{iComb};
            
            if iscell(EventsTestTable.(tempEventsTestFilterForColumn))
                matchIndicator = matchIndicator | ( strcmp(EventsTestTable.(tempEventsTestFilterForColumn), tempCompTest) );
            else
                matchIndicator = matchIndicator | ( EventsTestTable.(tempEventsTestFilterForColumn) ==  tempCompTest);
            end
        end
    end
    
    EventsTestTable = EventsTestTable(matchIndicator,:);
    
end
if nargout > 4
    EventsTestTableExcluded = EventsTestTable(~matchIndicator,:);
end

matchIndicator = ones(size(EventsTargetTable,1),1);
if ~isempty(cfg.EventsTargetFilterForColumns)
    matchIndicator = zeros(size(EventsTargetTable,1),1);
    for iCol = 1:numel(cfg.EventsTargetFilterForColumns)
        tempEventsTargetFilterForColumn = cfg.EventsTargetFilterForColumns{iCol};
        for iComb = 1:numel(cfg.EventsTargetFilterValues)
            %iComp = 1
            %tempCompTest = cfg.EventsTestFilterValues{iComb};
            tempCompTarget = cfg.EventsTargetFilterValues{iComb};
            
            if iscell(EventsTargetTable.(tempEventsTargetFilterForColumn))
                matchIndicator = matchIndicator | ( strcmp(EventsTargetTable.(tempEventsTargetFilterForColumn), tempCompTarget) );
            else
                matchIndicator = matchIndicator | ( EventsTargetTable.(tempEventsTargetFilterForColumn) ==  tempCompTarget);
            end
        end
    end
    
    EventsTargetTable = EventsTargetTable(matchIndicator,:);
end

if nargout > 5
    EventsTargetTableExcluded = EventsTargetTable(~matchIndicator,:);
end

nEventsTest = size(EventsTestTable,1);
nEventsTarget = size(EventsTargetTable,1);

groupByMapOverlap = containers.Map();
groupByMapNonOverlap = containers.Map();

groupByMapAllTest = containers.Map();
groupByMapAllTarget = containers.Map();
% groupByMapNonOverlapTarget = containers.Map();


if strcmp(cfg.MismatchIdenticalEvents,'yes') || strcmp(cfg.MismatchDuplicateTestTargetOrTargetTestMatchingEvents,'yes')
    
    IDmergeString = '#';
    temp_test_id = {''};
    temp_target_id = {''};
    for iComb = 1:numel(cfg.EventsTestIDColumns)
        tempCompTest = cfg.EventsTestIDColumns{iComb};
        tempCompTarget = cfg.EventsTargetIDColumns{iComb};
        
        %dsEventsTest.(tempCompTest)
        %dsEventsTarget.(tempCompTarget)
        
        if iscell(EventsTestTable.(tempCompTest))
            temp_test_id = strcat(temp_test_id,{IDmergeString},EventsTestTable.(tempCompTest));
        else
            temp_test_id = strcat(temp_test_id,{IDmergeString},num2str(EventsTestTable.(tempCompTest),'%-g'));
        end
        
        if iscell(EventsTargetTable.(tempCompTarget))
            temp_target_id = strcat(temp_target_id,{IDmergeString},EventsTargetTable.(tempCompTarget));
        else
            temp_target_id = strcat(temp_target_id,{IDmergeString},num2str(EventsTargetTable.(tempCompTarget),'%-g'));
        end
    end
    
    EventsTestTable.duplication_id = temp_test_id;
    EventsTargetTable.duplication_id = temp_target_id;
    
    temp_test_id = [];
    temp_target_id = [];
    
    duplicateIDMapTest = containers.Map();
    duplicateIDMapTarget = containers.Map();
    
end

overlap_test_target = [];
nonoverlap_test = [];



temp_overlap_collector_iterator = 1;
temp_overlap_collector = {};
temp_overlap_collector{temp_overlap_collector_iterator} = [];

temp_nonoverlap_collector_iterator = 1;
temp_nonoverlap_collector = {};
temp_nonoverlap_collector{temp_nonoverlap_collector_iterator} = [];

ft_progress('init', 'text', ['Please wait...']);

column_prefix_test = 'te_';
column_prefix_target = 'ta_';


columnNamesTestNew = strcat(column_prefix_test,EventsTestTable.Properties.VariableNames);
columnNamesTargetNew = strcat(column_prefix_target,EventsTargetTable.Properties.VariableNames);

% cut the column names because length of columns cannot exced 64 charecters
% with matlab TABLE
for iCol = 1:numel(columnNamesTestNew)
    if length(columnNamesTestNew{iCol}) > namelengthmax
        columnNamesTestNew{iCol} = [columnNamesTestNew{iCol}(1:(min(end,namelengthmax)-2)) '_X'];
    end
end
for iCol = 1:numel(columnNamesTargetNew)
    if length(columnNamesTargetNew{iCol}) > namelengthmax
        columnNamesTargetNew{iCol} = [columnNamesTargetNew{iCol}(1:(min(end,namelengthmax)-2)) '_X'];
    end
end

timeStartPar = toc(ttic);
progress_count = 0;
mismatchIndicator_target = ones(nEventsTarget,1);
for iEvTest = 1:nEventsTest
    %iEvTest = 1
    
    progress_count = progress_count + 1;
    tempTimerNow = toc(ttic);
    time_left_min = (nEventsTest-progress_count)*(((tempTimerNow - timeStartPar)/60)/progress_count);
    ft_progress(progress_count/nEventsTest, ['Processing test event %d (%d percent, matchchunk %d, mismatchchunk %d) of %d against %d targets, remains: ~%d:%02d min'], progress_count, fix(100*progress_count/nEventsTest), temp_overlap_collector_iterator, temp_nonoverlap_collector_iterator, nEventsTest, nEventsTarget, fix(time_left_min),fix((time_left_min-fix(time_left_min))*60) );  % show string, x=i/N
    
    
    eventTest = EventsTestTable(iEvTest,:);
    matchIndicator = ones(nEventsTarget,1);
    
    
    if strcmp(cfg.MismatchIdenticalEvents,'yes')
        matchIndicator = matchIndicator & ~(strcmp(eventTest.duplication_id,EventsTargetTable.duplication_id));
    end
    
    
    for iComb = 1:numel(cfg.EventsTestCompareColumns)
        %iComp = 1
        tempCompTest = cfg.EventsTestCompareColumns{iComb};
        tempCompTarget = cfg.EventsTargetCompareColumns{iComb};
        
        if iscell(eventTest.(tempCompTest))
            matchIndicator = matchIndicator & ( strcmp(eventTest.(tempCompTest), EventsTargetTable.(tempCompTarget)) );
        else
            matchIndicator = matchIndicator & ( eventTest.(tempCompTest) == EventsTargetTable.(tempCompTarget) );
        end
    end
    

    
    
    groupBy = 'group';
    for iGroup = 1:numel(cfg.EventsTestGroupSummaryByColumns)
        %iComp = 1
        tempCompTest = cfg.EventsTestGroupSummaryByColumns{iGroup};
        
        if iscell(eventTest.(tempCompTest))
            groupBy = [groupBy GroupByConcatString eventTest.(tempCompTest){:}];
        else
            groupBy = [groupBy GroupByConcatString num2str(eventTest.(tempCompTest))];
        end
    end
    
    if ~isKey(groupByMapOverlap,{groupBy})
        groupByMapOverlap(groupBy) = 0;
    end
    
    if ~isKey(groupByMapNonOverlap,{groupBy})
        groupByMapNonOverlap(groupBy) = 0;
    end
    
%     if ~isKey(groupByMapNonOverlapTarget,{groupBy})
%         groupByMapNonOverlapTarget(groupBy) = 0;
%     end
    
    if ~isKey(groupByMapAllTest,{groupBy})
        groupByMapAllTest(groupBy) = 0;
    end
    
    
    groupByMapAllTarget(groupBy) = sum(matchIndicator);
    %groupByMapNonOverlapTarget(groupBy) = sum(mismatchIndicator_target & ~matchIndicator);

    matchIndicatorTarget = matchIndicator;
    
    groupByMapAllTest(groupBy) = groupByMapAllTest(groupBy) + 1;
    
    if strcmp(cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow,'yes')
        matchIndicator = matchIndicator & ( (eventTest.(cfg.EventsTestTimePointColumn) + cfg.EventTestTimePointOffsetTime ) >= (EventsTargetTable.(cfg.EventsTargetTimePointColumn) + cfg.EventTargetTimeWindowOffsetTime) );
        matchIndicator = matchIndicator & ( (eventTest.(cfg.EventsTestTimePointColumn) + cfg.EventTestTimePointOffsetTime ) <= (EventsTargetTable.(cfg.EventsTargetTimePointColumn2) + cfg.EventTargetTimeWindowOffsetTime2) );
    else
        matchIndicator = matchIndicator & ( (eventTest.(cfg.EventsTestTimePointColumn) + cfg.EventTestTimePointOffsetTime ) <= ((EventsTargetTable.(cfg.EventsTargetTimePointColumn) + cfg.EventTargetTimeWindowOffsetTime) + cfg.EventTargetTimeWindowPostOffsetTime) );
        matchIndicator = matchIndicator & ( (eventTest.(cfg.EventsTestTimePointColumn) + cfg.EventTestTimePointOffsetTime ) >= ((EventsTargetTable.(cfg.EventsTargetTimePointColumn) + cfg.EventTargetTimeWindowOffsetTime) - cfg.EventTargetTimeWindowPreOffsetTime) );
    end
    %matchIndicator(2:3) = 1;
    
    
    eventTest.Properties.VariableNames = columnNamesTestNew;
    
    mismatchIndicator_target = mismatchIndicator_target & ~matchIndicator;

    
    if any(matchIndicator)
        tempDsEventsTarget = EventsTargetTable(matchIndicator,:);
        if strcmp(cfg.MismatchDuplicateTestTargetOrTargetTestMatchingEvents,'yes')
            
            temp_test_target_id_strings = strcat(eventTest.([column_prefix_test 'duplication_id']),EventsTargetTable.duplication_id(matchIndicator));
            temp_target_test_id_strings = strcat(EventsTargetTable.duplication_id(matchIndicator),eventTest.([column_prefix_test 'duplication_id']));
            
            temp_already_contained_test_first = isKey(duplicateIDMapTest,temp_test_target_id_strings);
            temp_already_contained_target_first = isKey(duplicateIDMapTarget,temp_target_test_id_strings);
            
            temp_add_to_overlapp_notcontainted_before_index = ~( temp_already_contained_test_first | temp_already_contained_target_first );
            
            if any(~temp_already_contained_test_first)
                temp_add_to_map = temp_test_target_id_strings(~temp_already_contained_test_first);
                duplicateIDMapTest = [duplicateIDMapTest;containers.Map(temp_add_to_map,ones(numel(temp_add_to_map),1))];
            end
            if any(~temp_already_contained_target_first)
                temp_add_to_map = temp_target_test_id_strings(~temp_already_contained_target_first);
                duplicateIDMapTarget = [duplicateIDMapTarget;containers.Map(temp_add_to_map,ones(numel(temp_add_to_map),1))];
            end
            tempDsEventsTarget = tempDsEventsTarget(temp_add_to_overlapp_notcontainted_before_index,:);
        end
        
        tempDsEventsTarget.Properties.VariableNames = columnNamesTargetNew;
        nOverlaps = size(tempDsEventsTarget,1);
        
        if nargout > 1
            for iTarRow = 1:nOverlaps
                
                if size(temp_overlap_collector{temp_overlap_collector_iterator},1) > 0
                    %if size(overlap,1) > 0
                    temp_overlap_collector{temp_overlap_collector_iterator}  = cat(1,temp_overlap_collector{temp_overlap_collector_iterator} ,cat(2,eventTest,tempDsEventsTarget(iTarRow,:)));
                    %overlap = cat(1,overlap,cat(2,eventTest,tempDsEventsTarget(iTarRow,:)));
                else
                    temp_overlap_collector{temp_overlap_collector_iterator} = cat(2,eventTest,tempDsEventsTarget(iTarRow,:));
                    %overlap = cat(2,eventTest,tempDsEventsTarget(iTarRow,:));
                end
                if size(temp_overlap_collector{temp_overlap_collector_iterator},1) > cfg.concatbuffersize
                    temp_overlap_collector_iterator = temp_overlap_collector_iterator + 1;
                    temp_overlap_collector{temp_overlap_collector_iterator} = [];
                end
                
            end
        end
        groupByMapOverlap(groupBy) = groupByMapOverlap(groupBy) + nOverlaps;
    else
        
        if nargout > 2
            if size(temp_nonoverlap_collector{temp_nonoverlap_collector_iterator},1) > 0
                %if size(nonoverlap,1) > 0
                temp_nonoverlap_collector{temp_nonoverlap_collector_iterator}  = cat(1,temp_nonoverlap_collector{temp_nonoverlap_collector_iterator} ,eventTest);
                %nonoverlap = cat(1,nonoverlap,eventTest);
            else
                temp_nonoverlap_collector{temp_nonoverlap_collector_iterator} = eventTest;
                %nonoverlap = eventTest;
            end
            if size(temp_nonoverlap_collector{temp_nonoverlap_collector_iterator},1) > cfg.concatbuffersize
                temp_nonoverlap_collector_iterator = temp_nonoverlap_collector_iterator + 1;
                temp_nonoverlap_collector{temp_nonoverlap_collector_iterator} = [];
            end
        end
        groupByMapNonOverlap(groupBy) = groupByMapNonOverlap(groupBy) + 1;
    end
end


if nargout > 1
    for iTemp_overlap_collector_iterator = 1:numel(temp_overlap_collector)
        if iTemp_overlap_collector_iterator == 1
            overlap_test_target = temp_overlap_collector{iTemp_overlap_collector_iterator};
        else
            overlap_test_target = cat(1,overlap_test_target,temp_overlap_collector{iTemp_overlap_collector_iterator});
        end
    end
    temp_overlap_collector = [];
end

if nargout > 2
    for iTemp_nonoverlap_collector_iterator = 1:numel(temp_nonoverlap_collector)
        if iTemp_nonoverlap_collector_iterator == 1
            nonoverlap_test = temp_nonoverlap_collector{iTemp_nonoverlap_collector_iterator};
        else
            nonoverlap_test = cat(1,nonoverlap_test,temp_nonoverlap_collector{iTemp_nonoverlap_collector_iterator});
        end
    end
    temp_nonoverlap_collector = [];
end

if nargout > 3
    nonoverlap_target = EventsTargetTable(mismatchIndicator_target,:);
end


%     if isempty(nonoverlap)
%         for iDScol = 1:numel(columnNamesTestNew)
%             nonoverlap = cat(2,nonoverlap,table([],'VariableNames',{num2str(iDScol)}));
%         end
%         nonoverlap = set(nonoverlap,'VariableNames',columnNamesTestNew);
%     end
%
%     if isempty(overlap)
%         for iDScol = 1:numel([columnNamesTestNew,columnNamesTargetNew])
%             overlap = cat(2,overlap,table([],'VariableNames',{num2str(iDScol)}));
%         end
%         overlap = set(overlap,'VariableNames',[columnNamesTestNew,columnNamesTargetNew]);
%     end

ft_progress('close');
varNames = {'used_2nd_column_and_offset_not_pre_and_post_and_offsets',[column_prefix_test 'compare_columns'], [column_prefix_target 'compare_columns'],[column_prefix_test 'timepoint_column'], [column_prefix_target 'timepoint_column'], [column_prefix_target 'offset'], [column_prefix_target 'pre_offset'], [column_prefix_target 'post_offset'], [column_prefix_test 'offset']};

tempEventsTestCompareColumns = cfg.EventsTestCompareColumns;
tempEventsTargetCompareColumns = cfg.EventsTargetCompareColumns;

if size(tempEventsTestCompareColumns,2) > 1
    tempEventsTestCompareColumns = {strjoin(cfg.EventsTestCompareColumns,' ')};
    tempEventsTargetCompareColumns = {strjoin(cfg.EventsTargetCompareColumns,' ')};
end


if nargout > 1
    %     if true
    addRightOverlap = table(repmat({cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow},size(overlap_test_target,1),1),...
        repmat(tempEventsTestCompareColumns,size(overlap_test_target,1),1),...
        repmat(tempEventsTargetCompareColumns,size(overlap_test_target,1),1),...
        repmat({cfg.EventsTestTimePointColumn},size(overlap_test_target,1),1),...
        repmat({cfg.EventsTargetTimePointColumn},size(overlap_test_target,1),1),...
        repmat(cfg.EventTargetTimeWindowOffsetTime,size(overlap_test_target,1),1),...
        repmat(cfg.EventTargetTimeWindowPreOffsetTime,size(overlap_test_target,1),1),...
        repmat(cfg.EventTargetTimeWindowPostOffsetTime,size(overlap_test_target,1),1),...
        repmat(cfg.EventTestTimePointOffsetTime,size(overlap_test_target,1),1),...
        'VariableNames',varNames);
    %     else
    %         addRightOverlap = table(repmat({cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow},size(overlap,1),1),repmat(tempEventsTestCompareColumns,size(overlap,1),1),...
    %             repmat(tempEventsTargetCompareColumns,size(overlap,1),1),...
    %             repmat({cfg.EventsTestTimePointColumn},size(overlap,1),1),repmat({cfg.EventsTargetTimePointColumn},size(overlap,1),1),...
    %             {repmat(cfg.EventTargetTimeWindowOffsetTime,size(overlap,1),1)},{repmat(cfg.EventTargetTimeWindowPreOffsetTime,size(overlap,1),1)},{repmat(cfg.EventTargetTimeWindowPostOffsetTime,size(overlap,1),1)},{repmat(cfg.EventTestTimePointOffsetTime,size(overlap,1),1)},'VariableNames',varNames);
    %     end
    
    overlap_test_target = cat(2,overlap_test_target,addRightOverlap);
    overlap_test_target = cat(2,table(repmat({eventsTestIDname},size(overlap_test_target,1),1),repmat({eventsTargetIDname},size(overlap_test_target,1),1),'VariableNames',{[column_prefix_test 'event_ori'],[column_prefix_target 'event_ori']}),overlap_test_target);

    overlap_test_target.test_minus_target_delay = overlap_test_target.([column_prefix_test cfg.EventsTestTimePointColumn]) - overlap_test_target.([column_prefix_target cfg.EventsTargetTimePointColumn]);

end

if nargout > 2
    %     if size(nonoverlap,1) > 1
    ncols = size(nonoverlap_test,1);
    addRightNonOverlap = table(repmat({cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow},ncols,1),...
        repmat(tempEventsTestCompareColumns,ncols,1),...
        repmat(tempEventsTargetCompareColumns,ncols,1),...
        repmat({cfg.EventsTestTimePointColumn},ncols,1),...
        repmat({cfg.EventsTargetTimePointColumn},ncols,1),...
        repmat(cfg.EventTargetTimeWindowOffsetTime,ncols,1),...
        repmat(cfg.EventTargetTimeWindowPreOffsetTime,ncols,1),...
        repmat(cfg.EventTargetTimeWindowPostOffsetTime,ncols,1),...
        repmat(cfg.EventTestTimePointOffsetTime,ncols,1),...
        'VariableNames',varNames);
    %     else
    %         addRightNonOverlap = table(repmat({cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow},size(nonoverlap,1),1),repmat(tempEventsTestCompareColumns,size(nonoverlap,1),1),...
    %             repmat(tempEventsTargetCompareColumns,size(nonoverlap,1),1),...
    %             repmat(cfg.EventsTestTimePointColumn,size(nonoverlap,1),1),repmat(cfg.EventsTargetTimePointColumn,size(nonoverlap,1),1),...
    %             {repmat(cfg.EventTargetTimeWindowOffsetTime,size(nonoverlap,1),1)},{repmat(cfg.EventTargetTimeWindowPreOffsetTime,size(nonoverlap,1),1)},{repmat(cfg.EventTargetTimeWindowPostOffsetTime,size(nonoverlap,1),1)},{repmat(cfg.EventTestTimePointOffsetTime,size(nonoverlap,1),1)},'VariableNames',varNames);
    %     end
    
    nonoverlap_test = cat(2,nonoverlap_test,addRightNonOverlap);
    nonoverlap_test = cat(2,table(repmat({eventsTestIDname},ncols,1),repmat({eventsTargetIDname},ncols,1),'VariableNames',{[column_prefix_test 'event_ori'],[column_prefix_target 'event_ori']}),nonoverlap_test);
end

if nargout > 3
    
    %     if size(nonoverlap,1) > 1
    ncols = size(nonoverlap_target,1);
    addRightNonOverlapTarget = table(repmat({cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow},ncols,1),...
        repmat(tempEventsTestCompareColumns,ncols,1),...
        repmat(tempEventsTargetCompareColumns,ncols,1),...
        repmat({cfg.EventsTestTimePointColumn},ncols,1),...
        repmat({cfg.EventsTargetTimePointColumn},ncols,1),...
        repmat(cfg.EventTargetTimeWindowOffsetTime,ncols,1),...
        repmat(cfg.EventTargetTimeWindowPreOffsetTime,ncols,1),...
        repmat(cfg.EventTargetTimeWindowPostOffsetTime,ncols,1),...
        repmat(cfg.EventTestTimePointOffsetTime,ncols,1),...
        'VariableNames',varNames);
    %     else
    %         addRightNonOverlap = table(repmat({cfg.UseSecondColumnAndOnlyOffsetsForTimeWindow},size(nonoverlap,1),1),repmat(tempEventsTestCompareColumns,size(nonoverlap,1),1),...
    %             repmat(tempEventsTargetCompareColumns,size(nonoverlap,1),1),...
    %             repmat(cfg.EventsTestTimePointColumn,size(nonoverlap,1),1),repmat(cfg.EventsTargetTimePointColumn,size(nonoverlap,1),1),...
    %             {repmat(cfg.EventTargetTimeWindowOffsetTime,size(nonoverlap,1),1)},{repmat(cfg.EventTargetTimeWindowPreOffsetTime,size(nonoverlap,1),1)},{repmat(cfg.EventTargetTimeWindowPostOffsetTime,size(nonoverlap,1),1)},{repmat(cfg.EventTestTimePointOffsetTime,size(nonoverlap,1),1)},'VariableNames',varNames);
    %     end
    
    nonoverlap_target = cat(2,nonoverlap_target,addRightNonOverlapTarget);
    nonoverlap_target = cat(2,table(repmat({eventsTestIDname},size(nonoverlap_target,1),1),repmat({eventsTargetIDname},ncols,1),'VariableNames',{[column_prefix_test 'event_ori'],[column_prefix_target 'event_ori']}),nonoverlap_target);
    
end


tempEventsTestGroupSummaryByColumns = cfg.EventsTestGroupSummaryByColumns;

if size(tempEventsTestGroupSummaryByColumns,2) > 1
    tempEventsTestGroupSummaryByColumns = {strjoin(cfg.EventsTestGroupSummaryByColumns,' ')};
end

groups = keys(groupByMapOverlap);
if size(groups,2) < 2
    summary = table(repmat('group',2,1),'VariableNames',{'group'});
    summary(2,:) = [];
else
    summary = table(repmat('group',size(groups,2),1),'VariableNames',{'group'});
end

for iGroup = 1:numel(cfg.EventsTestGroupSummaryByColumns)
    tempGroupTest = cfg.EventsTestGroupSummaryByColumns{iGroup};
    summary = cat(2,summary,table(repmat(cellstr('group'),size(groups,2),1),'VariableNames',{tempGroupTest}));
end

summary = cat(2,summary,table(cell2mat(values(groupByMapOverlap))','VariableNames',{[column_prefix_test 'match_grouped']}));
summary = cat(2,summary,table(cell2mat(values(groupByMapNonOverlap))','VariableNames',{[column_prefix_test 'mismatch_grouped']}));
% summary = cat(2,summary,table(cell2mat(values(groupByMapNonOverlapTarget))','VariableNames',{[column_prefix_target 'mismatch_ungrouped']}));
summary = cat(2,summary,table(cell2mat(values(groupByMapAllTest))','VariableNames',{[column_prefix_test 'grouped']}));
summary = cat(2,summary,table(cell2mat(values(groupByMapAllTarget))','VariableNames',{[column_prefix_target 'ungrouped']}));
summary = cat(2,summary,table(((summary.([column_prefix_test 'match_grouped']) + summary.([column_prefix_test 'mismatch_grouped'])) - summary.([column_prefix_test 'grouped'])),'VariableNames',{[column_prefix_test 'matches_' column_prefix_target 'more_than_once']}));
summary = cat(2,summary,table(repmat(tempEventsTestCompareColumns,size(summary,1),1),'VariableNames',{[column_prefix_test 'compare_columns']}));
summary = cat(2,summary,table(repmat(tempEventsTargetCompareColumns,size(summary,1),1),'VariableNames',{[column_prefix_target 'compare_columns']}));
summary = cat(2,summary,table(repmat(tempEventsTestGroupSummaryByColumns,size(summary,1),1),'VariableNames',{[column_prefix_test 'group_by_columns']}));
EventsTestFilterForColumns_temp = cfg.EventsTestFilterForColumns;
if isempty(cfg.EventsTestFilterForColumns)
    EventsTestFilterForColumns_temp = {''};
end
summary = cat(2,summary,table(repmat(EventsTestFilterForColumns_temp,size(summary,1),1),'VariableNames',{[column_prefix_test 'filter_column']}));
EventsTargetFilterForColumns_temp = cfg.EventsTargetFilterForColumns;
if isempty(cfg.EventsTargetFilterForColumns)
    EventsTargetFilterForColumns_temp = {''};
end

summary = cat(2,summary,table(repmat(EventsTargetFilterForColumns_temp,size(summary,1),1),'VariableNames',{[column_prefix_target 'filter_column']}));
%summary = cat(2,summary,table(repmat(cfg.EventsTestFilterValues,size(summary,1),1),'VariableNames',{[column_prefix_test 'filter_value']}));
%summary = cat(2,summary,table(repmat(cfg.EventsTargetFilterValues,size(summary,1),1),'VariableNames',{[column_prefix_target 'filter_value']}));

for iGroupKey = 1:size(groups,2)
    groupsSplits = strsplit(groups{iGroupKey},GroupByConcatString);
    for iGroup = 1:numel(cfg.EventsTestGroupSummaryByColumns)
        tempGroupTest = cfg.EventsTestGroupSummaryByColumns{iGroup};
        val = groupsSplits(iGroup+1);
        if ~iscell(res_test.table.(tempGroupTest))
            val = str2num(val{:});
            if iscell(summary.(tempGroupTest))
                summary.(tempGroupTest) = repmat(NaN,length(summary.(tempGroupTest)),1);
            end
        end
        summary.(tempGroupTest)(iGroupKey) = val;
    end
end

if numel(cfg.EventsTestGroupSummaryByColumns) > 0
    summary.group = [];
end



res_summary = [];
res_summary.ori = functionname;
res_summary.type = 'summary';
res_summary.cfg = cfg;
res_summary.table = summary;

if nargout > 1
    res_match_test_target = [];
    res_match_test_target.ori = functionname;
    res_match_test_target.type = 'match_test_target';
    res_match_test_target.cfg = cfg;
    res_match_test_target.table = overlap_test_target;
end

if nargout > 2
    res_mismatch_test = [];
    res_mismatch_test.ori = functionname;
    res_mismatch_test.type = 'mismatch_test';
    res_mismatch_test.cfg = cfg;
    if isempty(nonoverlap_test)
        colnames = cat(2,res_test.Properties.VariableNames,varNames);
        nonoverlap_test = cell2table(cell(0,numel(colnames)),'VariableNames',colnames);
    end
    res_mismatch_test.table = nonoverlap_test;
end

if nargout > 3
    res_mismatch_target = [];
    res_mismatch_target.ori = functionname;
    res_mismatch_target.type = 'mismatch_target';
    res_mismatch_target.cfg = cfg;
    if isempty(nonoverlap_target)
        colnames = cat(2,res_target.Properties.VariableNames,varNames);
        nonoverlap_target = cell2table(cell(0,numel(colnames)),'VariableNames',colnames);
    end
    res_mismatch_target.table = nonoverlap_target;
end

if nargout > 4
    res_excluded_test = [];
    res_excluded_test.ori = functionname;
    res_excluded_test.type = 'excluded_target';
    res_excluded_test.cfg = cfg;
    if isempty(EventsTestTableExcluded)
        colnames = res_test.table.Properties.VariableNames;
        EventsTestTableExcluded = cell2table(cell(0,numel(colnames)),'VariableNames',colnames);
    end
    res_excluded_test.table = EventsTestTableExcluded;
end

if nargout > 5
    res_excluded_target = [];
    res_excluded_target.ori = functionname;
    res_excluded_target.type = 'excluded_target';
    res_excluded_target.cfg = cfg;
    if isempty(EventsTargetTableExcluded)
        colnames = res_target.table.Properties.VariableNames;
        EventsTargetTableExcluded = cell2table(cell(0,numel(colnames)),'VariableNames',colnames);
    end
    res_excluded_target.table = EventsTargetTableExcluded;
end

fprintf([functionname ' function finished\n']);
toc
memtoc
end
