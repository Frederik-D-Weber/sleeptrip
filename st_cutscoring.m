function [scorings] = st_cutscoring(cfg, scoring)

% ST_CUTSCORING cuts a scoring into overlapping or non overlaping pieces
% of a specific length relative to sleeponset or other references or just
% cut out one specific part
%
% Use as
%   [scorings] = st_cutscoring(cfg,scoring)
%
%   scoring is a structure provided by ST_READ_SCORING
%   it returns many scoring sturctures stored in a cell array
%
% Configuration parameters are
%   cfg.length    = single number (in epochs)
%                   Note that the tail of the scoring that does not fit the
%                   requested length is discarded
%
% Alternatively on can define cfg.start and cfg.end epoch (inclusive) either in
%                       addition to cfg.length. 
%                       if cf.length is not provided it is inferred by
%                       cfg.start and cfg.end
%   cfg.start         = either single number or a Nx1 vector at which epoch(s) to start (in
%                       epochs) OR 'sleeponset' or 'lightsoff' (default = 'sleeponset')
%   cfg.end           = single number number or a Nx1 vector at which epoch(s) to end (in epochs)
%                       use Inf if you want to use always until the end of the
%                       available epochs (default = Inf)
%
%
% Configuration parameters are
%   cfg.overlap       = single number (in epochs) specifying the number of
%                       epochs that overlap, negative overlap are gaps (e.g. 0 = no overlap, -1 is one epoch gap) 
%                       (default = 0)
%   cfg.offset        = single number with the offset relative to the start
%                       (in epochs) (default = 0)
%   cfg.sleeponsetdef = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                       'NR' or 'NR' or 'XR', see ST_SCORINGDESCRIPTIVES for details (default = 'N1_XR')
%
% See also ST_READ_SCORING

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

% set the defaults
cfg.overlap                 = ft_getopt(cfg, 'overlap', 0);
cfg.start                   = ft_getopt(cfg, 'start', 'sleeponset');
cfg.offset                  = ft_getopt(cfg, 'offset', 0);
cfg.end                     = ft_getopt(cfg, 'end', Inf);
cfg.sleeponsetdef           = ft_getopt(cfg, 'sleeponsetdef', 'N1_XR');

if ~isfield(cfg,'length')
    cfg.length = cfg.end - cfg.start + 1;
end

if ~(cfg.overlap < cfg.length)
    ft_error('cfg.overlap = %d, but must be chosen smaller than cfg.length = %d',cfg.overlap,cfg.length)
end

hasLightsOff = false;
lightsOffMoment = 0;
if isfield(scoring, 'lightsoff')
    hasLightsOff = true;
    lightsOffMoment = scoring.lightsoff;
else
    lightsOffMoment = 0;
    ft_warning('The lights off moment was not provided in the scoring structure.\n The beginning of the scoring is thus assumed as lights off.');
end
startEpoch = 1;
hasMultipleStarts = false;
if ischar(cfg.start)
    switch cfg.start
        case 'sleeponset'
            [startEpochs preoffsetnumber onsetepoch] = st_sleeponset(cfg,scoring);
        case 'lightsoff'
            startEpochs = find(floor(lightsOffMoment/scoring.epochlength)+1,1);    
        otherwise
            ft_error('cfg.start option %s not known.',cfg.start)
    end
elseif isnumeric(cfg.start)
        if numel(cfg.start) > 1
            hasMultipleStarts = true;
        end
        if any(cfg.start < 1)
            ft_error('cfg.start must be an integer number greater or equal to 1')
        end
        if ~all(size(cfg.end) == size(cfg.start))
            ft_error('if cfg.start and cfg.end are numeric, they must be matching dimensions or of same size')
        end
        if any((cfg.end - cfg.start) < 1)
            ft_error('cfg.start must correspond to cfg.end and always be smaller')
        end
        startEpochs = cfg.start;
else
            ft_error('cfg.start option not known or it is not an integer number greater or equal to 1')
end

if ~(cfg.overlap < cfg.length)
    ft_error('cfg.overlap = %d, but must be chosen smaller than cfg.length = %d',cfg.overlap,cfg.length)
end


scorings = {};

endEpochs = cfg.end;

for iStart = 1:numel(startEpochs)
    
    snipLength = cfg.length(iStart);
    startEpoch = startEpochs(iStart);
    endEpoch = endEpochs(iStart);
    
    
    snipStart = startEpoch + cfg.offset;
    
    maxEpoch = min(numel(scoring.epochs),endEpoch);
    while (snipStart + snipLength - 1) <= maxEpoch
        snip = scoring;
        snipIndex = snipStart:(snipStart + snipLength - 1);
        snip.epochs = scoring.epochs(snipIndex);
        snip.excluded = scoring.excluded(snipIndex);
        if isfield(scoring,'numbers')
            snip.numbers = scoring.numbers(snipIndex);
        end
        if isfield(scoring,'prob')
            snip.prob = scoring.prob(:,snipIndex);
        end
        
        time_start = (snipStart-1) * scoring.epochlength + scoring.dataoffset;
        time_stop = max(snipIndex) * scoring.epochlength + scoring.dataoffset;
        if isfield(scoring, 'arousals')
            scoring.arousals = scoring.arousals((scoring.arousals.start >= time_start) & (scoring.arousals.start < time_stop),:);
        end
        if isfield(scoring, 'events')
            scoring.events = scoring.events((scoring.events.start >= time_start) & (scoring.events.start < time_stop),:);
        end

        snip.dataoffset = time_start;
        if hasLightsOff
            snip.lightsoff = lightsOffMoment - ((snipStart - 1) * scoring.epochlength);
        end
        scorings = cat(2,scorings,{snip});
        
        snipStart = snipStart + snipLength - cfg.overlap;
    end
end



