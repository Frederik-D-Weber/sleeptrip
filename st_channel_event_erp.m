function [timelock reschannelcolumnname channels] = st_channel_event_erp(cfg, res_event, data)

% ST_CHANNEL_EVENT_ERP creates the event related potential of the average 
% event for each channel separately form a res structure
% (e.g. ST_SPINDLES, ST_SLOWWAVES, ST_COOC). It also can determin the
% reschannelcolumnname that was used as an output, which is handy for
% subsequent functions to keep track of, as this might change from res to
%
% Use as
%   [timelock] = st_channel_event_erp(cfg, res, data)
%   [timelock reschannelcolumnname] = st_channel_event_erp(cfg, res, data)
%    if cfg.keeptrials = 'yes'
%   [timelocks reschannelcolumnname channels] = st_channel_event_erp(cfg, res_event, data)
%   [timelocks channels] = st_channel_event_erp(cfg, res_event, data)
%    where timelock is a cell with each event for each of the channels
%
% Required configuration parameters are:
%   cfg.eventtimecolumn = string, giving the column names to hold the even times = 0 in the in res.table (default = 'seconds_trough_max')
%   cfg.eventtimeoffset = number as an offset for defining where time = 0 is (default = 0)
%   cfg.bounds          = the number of seconds around the event time plus
%                         offset, e.g. +-5 seconds around the eventtime = 0 (default = [-5 5])
%   cfg.baselinecorrect = if to correct the signal by a baseline (see FT_TIMELOCKBASELINE for details)(default = 'yes');
%   cfg.baseline        = time of the baseline window with respect to time = 0 (see FT_TIMELOCKBASELINE for details) (default = cfg.bounds);
%   cfg.maxevents       = number of trials/event to maximally use per channel before sampled out (by random) (default is unset and unlimited trial number)
%   cfg.randseed        = numer as a seed for the randomization functions.
%   cfg.eventsamplemethod = the method to select the maximal number (N) of events in case cfg.maxevents is set, either 'first' (for the first N) 'random' or 'last' (the last N), (default = 'random').
%   cfg.keeptrials      = 'yes' or 'no', return individual trials or average (default = 'no')
%
% Optional configuration parameters from ST_RES_FILTER are:
%
%   cfg.channel       = the channels you want to select (default = 'all');
%   cfg.reschannelcolumnname = name of column in the res.table (default is
%                        chosen by the firs column that contains the string
%                        'channel';
%   cfg.filtercolumns = either a string or a cellstr of the columns you want to filter for, e.g.
%                        {'resnum', 'freq'}
%   cfg.filtervalues  = either a string/value or cell of cells/cellstr, with the corresponding values to
%                       the columns defined cfg.filtercolumns, e.g.
%                        {{1}, {4:5}}
%
%
% See also ST_SPINDLES, ST_SLOWWAVES, ST_RES_FILTER, FT_TIMELOCKBASELINE, FT_PREPROCESSING, FT_APPLY_MONTAGE

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

channels = {};
% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
hasMaxEvents = false;
if isfield(cfg, 'maxevents')
    hasMaxEvents = true;
end
cfg.bounds   = ft_getopt(cfg, 'bounds', [-5 5]);

cfg.baselinecorrect = ft_getopt(cfg, 'baselinecorrect', 'no');
cfg.baseline   = ft_getopt(cfg, 'baseline', cfg.bounds);
cfg.baselinetype   = ft_getopt(cfg, 'baselinetype', 'normchange');
cfg.eventtimecolumn   = ft_getopt(cfg, 'eventtimecolumn', 'seconds_trough_max');
cfg.eventtimeoffset   = ft_getopt(cfg, 'eventtimeoffset', 0);
cfg.keeptrials   = ft_getopt(cfg, 'keeptrials', 'no');


if istrue(cfg.baselinecorrect)
bounds = [min(min(cfg.baseline),cfg.bounds(1)) max(max(cfg.baseline),cfg.bounds(2))];
else
    bounds = cfg.bounds;
end

cfg.randseed   = ft_getopt(cfg, 'randseed', 42);
cfg.eventsamplemethod = ft_getopt(cfg, 'eventsamplemethod', 'random');

data = ft_checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

fprintf([functionname ' function initialized\n']);

[res_event reschannelcolumnname] = st_res_filter(cfg, res_event);

cfg_sd = [];
cfg_sd.channel = cfg.channel;
data = ft_selectdata(cfg_sd, data);

nChannels = numel(data.label);

event_time_chs = cell(1,nChannels);

for iCh = 1:nChannels
    
    ch = data.label{iCh};
    
    chanrows = strcmp(res_event.table.(reschannelcolumnname),{ch});
    eventsSeconds = res_event.table.(cfg.eventtimecolumn);
    eventsSeconds = eventsSeconds(chanrows) + cfg.eventtimeoffset;
    
    if hasMaxEvents
        switch cfg.eventsamplemethod
            case 'random'
                ind_sample = randsample(numel(eventsSeconds),min(cfg.maxevents,numel(eventsSeconds)),false);
                ind_sample = sort(ind_sample);
            case 'first'
                ind_sample = 1:min(numel(eventsSeconds),cfg.maxevents);
            case 'last'
                ind_sample = 1:min(numel(eventsSeconds),cfg.maxevents);
                ind_sample = (numel(eventsSeconds) - flip(ind_sample))+1;
        end
        eventsSeconds = eventsSeconds(ind_sample);
    end
    
    %check if there are any events
    if ~isempty(eventsSeconds)
        cfg_sd = [];
        cfg_sd.channel = ch;
        data_ch = ft_selectdata(cfg_sd, data);
        
        cfg_rt = [];
        cfg_rt.seconds = eventsSeconds;
        cfg_rt.bounds = bounds; % 10 seconds around the events
        data_events_ch = ft_redefinetrial(cfg_rt, data_ch);
        
        if ~isempty(data_events_ch.trial)
            %View the event average signal timelocked to the trough.
            cfg_tl        = cfg;
            cfg_tl.channel = ch;
            event_timelock_ch = ft_timelockanalysis(cfg_tl, data_events_ch);
                       
            if istrue(cfg.baselinecorrect)
                % view the time-frequency of the event, for each channel
                cfg_bl                = cfg;
                cfg_bl.baseline       = cfg.baseline;
                event_timelock_ch     = ft_timelockbaseline(cfg, event_timelock_ch);
            end
            
            cfg_sd = [];
            cfg_sd.latency = cfg.bounds;
            event_timelock_ch = ft_selectdata(cfg_sd,event_timelock_ch);
            
            event_time_chs{iCh} = event_timelock_ch;
            channels = cat(1,channels,ch);
        end
    else
        ft_warning('No events left for channel %s',ch)
    end
    
end

idx_empty_channels = cellfun(@isempty, event_time_chs);
event_time_chs = event_time_chs(~idx_empty_channels);

if all(idx_empty_channels)
    timelock = [];
    ft_warning('all channels produced had empty timelock structures.')
elseif any(idx_empty_channels)
    ft_warning('at least one channel produced had empty timelock structures and was thus excluded from the resutling timelock structure.')
else
    if istrue(cfg.keeptrials)
        timelock = event_time_chs;
    else
        timelock = ft_appendtimelock([], event_time_chs{:});
    end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance timelock
ft_postamble history    timelock
ft_postamble savevar    timelock


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
