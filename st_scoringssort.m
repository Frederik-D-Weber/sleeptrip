function [scorings_sorted sums_sorted sums_ori idx_ori] = st_scoringssort(cfg, scorings, varargin)

% ST_SCORINGSSORT sorts a cell structure with muliple 
%
% Use as
%   scorings_sorted                = st_scoringssort(cfg, scorings, res_cycles)
%   [scorings_sorted idx_ori]      = st_scoringssort(cfg, scorings, res_cycles)
%   [scorings_sorted sums sums_ori idx_ori] = st_scoringssort(cfg, scorings, res_cycles)
%   scorings_sorted                = st_scoringssort(cfg, scorings)
%   [scorings_sorted idx_ori]      = st_scoringssort(cfg, scorings)
%   [scorings_sorted sums sums_ori idx_ori] = st_scoringssort(cfg, scorings)
%
%
% Optional configuration parameters are:
%   cfg.sortby        = string with option to sort as follows
%              'cycles': chosen by default if res_cycles is
%                        provided as an input 
%                        see:
%                        cfg.cycles and 
%                        cfg.cycleproperty
%       'excludedcount': by count of excluded epochs in scoring
%          'epochcount': by count of all epochs in scoring
%  'epochcountunscored': by count of epochs unscored/unknown in scoring
%       'sleepunscored': by count of epochs unscored/unknown in scoring
%                        within sleep.
%       'sleepduration': by duration of sleep (i.e. total sleep time/duration)
%          'sleeponset': by the start of sleep onset from scoring beginning
%  'sleeponsetduration': by the start of sleep onset from lightsoff
%         'sleepoffset': by the start of sleep offset (after last sleep epoch) from scoring beginning
%'sleeponsettomidsleep': by the time from sleep onset to mid sleep time
%      'timetomidsleep': by the time from first epoch to mid sleep time
%    'midsleeptooffset': by the time from mid sleep to sleep offset (after
%                        last sleep epoch
%           'lightsoff': by the time of lights-off if provided
%        'descriptives': by st_scoringdescriptives output, 
%                        see cfg.descriptive for details
%                      (default = 'epochcount')
%   cfg.sortdir       = string, either 'ascend' or 'descend' (default = 'ascend')
%                       see 'sort' function documentation for more details
%   cfg.descriptive   = string, only used if cfg.sortby = 'descriptive' 
%                       then this provides the column name as in 
%                       the result output from res.table from
%                       ST_SCORINGDESCRIPTIVES with numeric values
%                       e.g. 'Sleep_Onset_min'
%                       or 'SWS_onset_min'
%                       or 'SWS_min'
%                       or 'R_onset_min'
%                       or 'R_min'
%                       or 'N2_onset_min'
%                       or 'N2_min'
%                       this can also be muliple columns, then describe
%                       them like e.g. {'R_min','SWS_min'} and the columns
%                       are summed together (here the sum of REM and SWS in
%                       minutes.
%                      (default = 'Total_sleep_time_min')
%   cfg.sleeponsetdef = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                         'NR' or 'N2R' or 'XR' or 'AASM' or 'X2R' or 
%                         'N2' or 'N3' or 'SWS' or 'S4' or 'R', 
%                         see below for details (default = 'N1_XR') see
%                         ST_SLEEPONSET for details
%   cfg.cycles        = value or vector of values defining which sleep cycles to sort for
%                       e.g. [1 3] to sort for the first and third cycle if present
%                      (default = 1) i.e. the first cycle.
%   cfg.cycleproperty = string, must be a column in the res_cycles{:}.table
%                       either 'startepoch' 'endepoch' 'NRstartepoch'
%                       'NRendepoch' 'Rstartepoch' 'Rendepoch' 
%                       'durationepochs' 'durationNRepochs' 'durationRepochs'
%                       (default = 'durationepochs')
%
%
% See also ST_READ_SCORING, ST_SLEEPCYCLES, ST_SLEEPONSET, ST_CUTSCORING, ST_SCORINGDESCRIPTIVES

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

% set defaults
cfg.cycleproperty =   ft_getopt(cfg, 'cycleproperty', 'durationepochs');                                  
cfg.cycles        =   ft_getopt(cfg, 'cycles', 1);
cfg.sortdir       =   ft_getopt(cfg, 'sortdir', 'ascend');
cfg.sortby        =   ft_getopt(cfg, 'sortby', 'epochcount');
cfg.descriptive   =   ft_getopt(cfg, 'descriptive', 'Total_sleep_time_min');



fprintf([functionname ' function initialized\n']);

if ~iscell(scorings)
    scorings = {scorings};
end

%check consistency
if numel(scorings) > 1
    for iScoring = 2:numel(scorings)
        
        if (scorings{1}.epochlength ~= scorings{iScoring}.epochlength)
            ft_error(['epochlength of scoring number ' num2str(iScoring) ' inconsistent and does not match the first one.']);
        end
        
        if ~strcmp(scorings{1}.standard,scorings{iScoring}.standard)
            ft_warning(['scoring standard of scoring number ' num2str(iScoring) ' inconsistent and does not match the first one.']);
        end
    end
end

if nargin > 2
    res_cycles = varargin{1};
    if ~iscell(res_cycles)
        res_cycles = {res_cycles};
    end
    if numel(scorings) ~= numel(res_cycles)
        ft_error('Number of scorings (%d) must agree with number of res_cycles (%d)',numel(scorings),numel(res_cycles))
    end
    if ~strcmp(cfg.sortby,'cycles')
        ft_warning('Even though res_cycles is provided as input, it is not used for sorting.\n if you want to use set cfg.sortby = ''cycles'' and maybe change the cfg.cycleproperty setting.');
    end
    cfg.sortby = 'cycles';
elseif nargin <= 2
    if strcmp(cfg.sortby,'cycles')
        ft_error('cannot sort by cycles if no res_cycles cell structure is provided as input.\nPlease provide a res_cycle structure like the output from ST_SLEEPCYCLES')
    end
    cfg.sortby = 'epochcount';
end

sums_ori = nan(1,numel(scorings));
for iScoring = 1:numel(scorings)
    scoring = scorings{iScoring};
    switch cfg.sortby
        case 'cycles'
            res_cycle = res_cycles{iScoring};
            sums_ori(iScoring) = nansum(res_cycle.table(ismember(cfg.cycles,res_cycle.table.cycle),:).(cfg.cycleproperty));
        case 'excludedcount'
            sums_ori(iScoring) = sum(scoring.excluded);
        case 'epochcount'
            sums_ori(iScoring) = numel(scoring.epochs);
        case 'epochcountunscored'
            unknown = '?';
            if strcmp(scoring.standard, 'number')
                unknown = '-1';   
            end
            if isfield(scoring,'cfg')
                if isfield(scoring.cfg,'scoremap')
                    unknown = scoring.cfg.scoremap.unknown;
                end
            end
            sums_ori(iScoring) = sum(~strcmp(scoring.epochs,unknown));
        case 'sleepunscored'
            [onsetnumber, lastsleepstagenumber] = st_sleeponset(cfg,scoring);
            unknown = '?';
            if strcmp(scoring.standard, 'number')
                unknown = '-1';   
            end
            if isfield(scoring,'cfg')
                if isfield(scoring.cfg,'scoremap')
                    unknown = scoring.cfg.scoremap.unknown;
                end
            end
            sums_ori(iScoring) = sum(~strcmp(scoring.epochs(onsetnumber:lastsleepstagenumber),unknown));
        case 'sleepduration'
            [onsetnumber, lastsleepstagenumber] = st_sleeponset(cfg,scoring);
            sums_ori(iScoring) = (lastsleepstagenumber-onsetnumber+1);
        case 'sleeponset'
            [onsetnumber, lastsleepstagenumber] = st_sleeponset(cfg,scoring);
            sums_ori(iScoring) = onsetnumber;
        case 'sleeponsetduration'
            [onsetnumber, lastsleepstagenumber] = st_sleeponset(cfg,scoring);
            lightsoff = 0;
            if isfield(scoring, 'lightsoff')
                lightsoff = scoring.lightsoff;
            else
                ft_warning('scoring.lightsoff for scoing number %d not found.\nWill set scoring.lightsoff = 0',iScoring);
                scoring.lightsoff = lightsoff;
                scorings{iScoring} = scoring;
            end
            sums_ori(iScoring) = scoring.epochlength*(onsetnumber-1) - lightsoff;
        case 'sleepoffset'
            [onsetnumber, lastsleepstagenumber] = st_sleeponset(cfg,scoring);
            sums_ori(iScoring) = lastsleepstagenumber;
        case 'sleeponsettomidsleep'
            [onsetnumber, lastsleepstagenumber] = st_sleeponset(cfg,scoring);
            sums_ori(iScoring) = (lastsleepstagenumber-onsetnumber+1)/2;
        case 'timetomidsleep'
            [onsetnumber, lastsleepstagenumber] = st_sleeponset(cfg,scoring);
            sums_ori(iScoring) = onsetnumber-1+(lastsleepstagenumber-onsetnumber+1)/2;
        case 'midsleeptooffset'
            [onsetnumber, lastsleepstagenumber] = st_sleeponset(cfg,scoring);
            sums_ori(iScoring) = lastsleepstagenumber + 1 - (onsetnumber-1+(lastsleepstagenumber-onsetnumber+1)/2);
        case 'lightsoff'
            lightsoff = 0;
            if isfield(scoring, 'lightsoff')
                lightsoff = scoring.lightsoff;
            else
                ft_warning('scoring.lightsoff for scoing number %d not found.\nWill set scoring.lightsoff = 0',iScoring);
                scoring.lightsoff = lightsoff;
                scorings{iScoring} = scoring;
            end
            sums_ori(iScoring) = lightsoff;
        case 'descriptives'
            res_desc = st_scoringdescriptives(cfg,scoring);
            value = res_desc.table(:,cfg.descriptive);
            sums_ori(iScoring) = sum(value);
        otherwise
            ft_error('cfg.sortby = %s is unknown, see the help for valid options.', cfg.sortby)
    end
end
[sums_sorted, idx_ori] = sort(sums_ori,cfg.sortdir);
scorings_sorted = scorings(idx_ori);


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
