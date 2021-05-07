function [scoring] = st_exclude_events_scoring(cfg, scoring, varargin)

% ST_EXCLUDE_EVENTS_SCORING excludes the epochs in a scoring that overlap
% or fall together with events.
%
% Use as
%   [scoring] = st_exclude_events_scoring(cfg, scoring, event_seconds)
%   [scoring] = st_exclude_events_scoring(cfg, scoring, event_starts_from_dataonset_seconds, event_stops_from_dataonset_seconds)
%
%
% The configuration structure can specify
%   cfg.timebuffer           = a 1x2 number vector with the time buffer
%                              left and right of the event or the event
%                              bounds, respectively (default = [0 0], e.g. [-3 3] to extend the time points by a fixed 3 seconds in each time direction)
%
% See also ST_READ_SCORING

% Copyright (C) 2020-, Frederik D. Weber
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar scoring
ft_preamble provenance scoring
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
    return
end

% set the defaults
cfg.timebuffer      = ft_getopt(cfg, 'timebuffer', [0 0]);


isBoundayEvent = false;
starts = varargin{1};

if nargin > 3
    isBoundayEvent = true;
    stops = varargin{2};
end


nEvents = numel(starts);
nEpochs = numel(scoring.excluded);


if nEvents > 0
starts = starts - scoring.dataoffset;
if isBoundayEvent
    stops = stops - scoring.dataoffset;
else
    stops = starts;
end


starts = starts + cfg.timebuffer(1);
stops = stops + cfg.timebuffer(2);


starts = floor(starts/scoring.epochlength)+1;
stops = floor(stops/scoring.epochlength)+1;

exind = [];
for iEvent = 1:nEvents
    exind = [exind starts(iEvent):stops(iEvent)];
end

exind = unique(exind);

exind((exind < 0) | (exind > nEpochs)) = [];
scoring.excluded(exind(:)') = true;

end

            
% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous scoring
ft_postamble provenance scoring
ft_postamble history    scoring
ft_postamble savevar    scoring
end
