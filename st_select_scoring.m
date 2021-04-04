function [cfg ends] = st_select_scoring(cfg, data)

% ST_SELECT_SCORING selects trials of sleep stages and returns them in a
% "trl" structure like FT_DEFINETRIAL
% this can be used to extract the parts of and conserve the time in the original data
%
% Use as
%   [data] = st_select_scoring(cfg, data)
%   [cfg]  = st_select_scoring(cfg)
%   [begins ends] = st_select_scoring(cfg, data)
%
%  with the necessary parameters for the configuration
%   cfg.scoring  = a scoring structure as defined in ST_READ_SCORING
%   cfg.stages   = either a string or a Nx1 cell-array with strings that indicate sleep stages
%                  if no data structure is provided as input next to configuration also
%                  provide this parameter
%
%  if no data structure is provided as parameter then you need to define
%   cfg.dataset  = string with the filename
% 
%
% See also ST_READ_SCORING, ST_PREPROCESSING, FT_DEFINETRIAL

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


%assure the stages are a cellstr structure, even with one element
if ~iscellstr(cfg.stages)
    if ischar(cfg.stages)
        cfg.stages = {cfg.stages};
    else
        ft_error('the stages of interest must be defines as a string or a cellstr.');
    end
end
  
hasdata = false;
if nargin > 1
    % data structure provided
    % check if the input data is valid for this function, the input data must be raw
    data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'hassampleinfo', 'yes');
    if isfield(data, 'trial') && isfield(data, 'time')
        % check if data structure is likely continous data
        if ((numel(data.trial) ~= 1) || ~all(size(data.sampleinfo) == [1 2])) && ~(nargout > 1)
            ft_error('data structure does not look like continous data and has more than one trial')
        end
    end
    hasdata = true;
    if iscell(data.time)
        nSamplesInData = numel(data.time{1});
    else
        nSamplesInData = numel(data.time);
    end
    fsample = data.fsample;
else
    % no data structure provided
    hdr = ft_read_header(cfg.dataset);
    nSamplesInData = hdr.nSamples;
    fsample = hdr.Fs;
end


sleepStagesOfInterst = cfg.stages;
epochs = cfg.scoring.epochs;
excluded = cfg.scoring.excluded;

epochLengthSamples = cfg.scoring.epochlength*fsample;
offsetSamples = cfg.scoring.dataoffset*fsample;

if nargout <= 1
    nPossibleEpochsInData = floor((nSamplesInData-offsetSamples)/epochLengthSamples);
    nEpochsInScoring = numel(epochs);
    if (nEpochsInScoring > nPossibleEpochsInData)
        epochs = epochs(1:nPossibleEpochsInData);
        excluded = excluded(1:nPossibleEpochsInData);
        ft_warning('only %d epochs possible in data with epoch length of %d and scoring offset of %f seconds \nbut scoring has %d epochs.\nPlease check if scoring matches to the data.\nScoring was shortened to fit data, and some epochs were discarded.',nPossibleEpochsInData,cfg.scoring.epochlength,cfg.scoring.dataoffset,nEpochsInScoring);
    end
end

hypnEpochs = 1:numel(epochs);
hypnEpochsBeginsSamples = (((hypnEpochs - 1) * epochLengthSamples) + 1)';
hypnEpochsEndsSamples = (hypnEpochs * epochLengthSamples)';

epochsOfInterst = hypnEpochs(ismember(epochs',sleepStagesOfInterst) & (excluded' == 0));

cfg.nEpochs = length(epochsOfInterst);

begins = [];
ends = [];
    
if ~isempty(epochsOfInterst)
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(epochsOfInterst,1);
    begins = hypnEpochsBeginsSamples(consecBegins);
    ends = hypnEpochsEndsSamples(consecEnds);
end


if hasdata
cfg.contbegsample = begins+offsetSamples;
cfg.contendsample = ends+offsetSamples;


% requested sleep stages are not present in dat
if isempty(cfg.contbegsample)
    ft_warning('requested sleep stages are not present in the scoring, will return empty trials.');
end

else
cfg.trl = [begins+offsetSamples ends+offsetSamples begins-1+offsetSamples];
% requested sleep stages are not present in dat
if isempty(cfg.trl)
    ft_warning('requested sleep stages are not present in the scoring, will result in empty trials.');
end

end
ends = [];
if hasdata
    if nargout > 1
        if ~isempty(cfg.contbegsample)
            %             ends = (cfg.trl(:,2)-1)/data.fsample;
            %             cfg = (cfg.trl(:,1)-1)/data.fsample; % cfg = begins
            %               ends = data.time{1}(cfg.trl(:,2));
            %               cfg = data.time{1}(cfg.trl(:,1)); % cfg = begins
            ends = data.time{1}(cfg.contendsample)';
            cfg = data.time{1}(cfg.contbegsample)'; % cfg = begins
        else
            cfg = [];
        end
    else
        if isempty(cfg.contbegsample)
            data.time = {};
            data.trial = {};
            data.sampleinfo = [];
        else
            data = ft_redefinetrial(cfg, data); 
        end
        %FIXME set data.cfg.prev
        cfg = data;
    end
end