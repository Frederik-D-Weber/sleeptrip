function [res_event reschannelcolumnname] = st_channel_event_phase(cfg, res_event, data)

% ST_CHANNEL_EVENT_PHASE creates the event related potential of the average 
% event for each channel separately form a res structure
% (e.g. ST_SPINDLES, ST_SLOWWAVES, ST_COOC). It also can determin the
% reschannelcolumnname that was used as an output, which is handy for
% subsequent functions to keep track of, as this might change from res to
%
% Use as
%   [res] = st_channel_event_phase(cfg, res, data)
%   [res reschannelcolumnname] = st_channel_event_phase(cfg, res, data)
%
% Optional configuration parameters are:
%   cfg.eventtimecolumn = string, giving the column names to hold the even times = 0 in the in res.table (default = 'seconds_trough_max')
%   cfg.grammer         = string with data operations, SEE ST_DATAGRAMMER for details (default = 'hp 0.5 lp 2' for slow oscillation phase)
%   cfg.phasetype       = the type of phase and value range as
%                         'phase_wrapped_rad', 'phase_wrapped_turns'
%                         'phase_wrapped_deg', 'phase_unwrapped_rad',
%                         'phase_unwrapped_turns', 'phase_unwrapped_deg' ,
%                         SEE ST_DATAGRAMMER for details (default = 'phase_wrapped_rad')
%
%   cfg.eventtimeoffset = number as an offset for defining where time = 0 is (default = 0)
%   cfg.bounds          = the number of seconds around the event time plus
%                         offset, e.g. +- 2 samples around the eventtime = 0 (default = plut or minus 2 samples)
%
% Optional configuration parameters from ST_RES_FILTER are:
%
%   cfg.channel       = the channels you want to select (default = 'all');
%   cfg.reschannelcolumnname = name of column in the res.table (default is
%                        chosen by the first column that contains the string
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

% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.baselinecorrect = 'no';
cfg.keeptrials = 'yes';
cfg.bounds   = ft_getopt(cfg, 'bounds', [-2/data.fsample 2/data.fsample]);
cfg.eventtimecolumn   = ft_getopt(cfg, 'eventtimecolumn', 'seconds_trough_max');
cfg.eventtimeoffset   = ft_getopt(cfg, 'eventtimeoffset', 0);

cfg.grammer         = ft_getopt(cfg, 'grammer', 'hp 0.5 lp 2');
cfg.phasetype         = ft_getopt(cfg, 'phasetype', 'phase_wrapped_rad');

if istrue(cfg.baselinecorrect)
bounds = [min(min(cfg.baseline),cfg.bounds(1)) max(max(cfg.baseline),cfg.bounds(2))];
else
    bounds = cfg.bounds;
end

%delete bordering samples from the result structure.
res_event.table((res_event.table.(cfg.eventtimecolumn) + bounds(2) + 1/data.fsample) >= max(data.time{1}),:) = [];
res_event.table((res_event.table.(cfg.eventtimecolumn) + bounds(1) - 1/data.fsample) <= min(data.time{1}),:) = [];

cfg.randseed   = ft_getopt(cfg, 'randseed', 42);
cfg.eventsamplemethod = ft_getopt(cfg, 'eventsamplemethod', 'random');

data = ft_checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

fprintf([functionname ' function initialized\n']);

[res_event reschannelcolumnname] = st_res_filter(cfg, res_event);

cfg_sd = [];
cfg_sd.channel = cfg.channel;
data = ft_selectdata(cfg_sd, data);



%%% one grammer to multiple or single channels
cfg_dg = [];
cfg_dg.grammer = [cfg.grammer ' ' cfg.phasetype];
cfg_dg.channel = cfg.channel;
data = st_datagrammer(cfg_dg, data);


[timelocks reschannelcolumnname channels] = st_channel_event_erp(cfg, res_event, data);

res_event_table_phase_column_name = [cfg.eventtimecolumn '_' cfg.phasetype];
res_event.table.(res_event_table_phase_column_name) = nan(size(res_event.table,1),1);

for iTimelock = 1:numel(timelocks)
    timelock = timelocks{iTimelock};
    time0ind = find(timelock.time == min(abs(timelock.time)));
    time0vals = timelock.trial(:,time0ind);
    ch = channels{iTimelock};
    res_event.table.(res_event_table_phase_column_name)(ismember(res_event.table.(reschannelcolumnname),{ch})) = time0vals;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance res_event
ft_postamble history    res_event
ft_postamble savevar    res_event


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
