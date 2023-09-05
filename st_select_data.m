function [data] = st_select_data(cfg, data)

% ST_SELECT_DATA selects data trials of requested sleep stages from a continuous recording,
% and returns them in a "trl" structure like FT_DEFINETRIAL. If the scoring has epochs labeled for exclusion (e.g. from artifact detection),
% these are excluded by default. Alternatively, select data using options for FT_REDEFINETRIAL.
%
% Use as
%   [data] = st_select_data(cfg, data)
%
%  When providing a scoring:
%   cfg.scoring  = a scoring structure as defined in ST_READ_SCORING
%
%  optional parameters:
%   cfg.stages   = stages to select. either a string or a Nx1 cell-array of strings (default: all available stages)
%   cfg.usescoringexclusion = 'yes' or 'no'. Whether or not to consider the scoring's exclusion information (default: 'yes')
%
%  When NOT providing a scoring, this function mostly acts like a wrapper around FT_REDEFINETRIAL. Provide (only one of) the following:
%   cfg.contbegtime = single number or Nx1 vector, expressed in seconds relative to the start of the data
%   cfg.contendtime = single number or Nx1 vector, expressed in seconds relative to the start of the data
%
%   cfg.length    = segment the data into equal duration trials. single number (in seconds) of the required snippets, e.g. 60
%   cfg.overlap   = only to combine with cfg.length: single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
%
%   cfg.contbegsample = single number or Nx1 vector, expressed in samples relative to the start of the data
%   cfg.contendsample = single number or Nx1 vector, expressed in samples relative to the start of the data
%
%   Additional selection options understood by FT_REDEFINETRIAL may also be passed.
%
%   GENERAL optional parameters:
%   cfg.minlength = minimal trial duration needed for selection (default: 30 seconds)
%
% See also ST_READ_SCORING, ST_PREPROCESSING, ST_SCORING_ARTIFACT_LEVEL, FT_REDEFINETRIAL

% Copyright (C) 2022-, Roy Cox, Frederik D. Weber
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

%---input checks and defaults----
%ft_checkconfig(cfg,'required',{'scoring'});
cfg.usescoringexclusion = ft_getopt(cfg, 'usescoringexclusion', 'yes');
cfg.minlength = ft_getopt(cfg, 'minlength',30); %default 30s (typical epoch length)

% ---check scoring and stages----
hasScoring = false;

if isfield(cfg,'scoring')
    scoring = cfg.scoring;
    hasScoring = true;

    %stages only makes sense in presence of scoring
    cfg.stages = ft_getopt(cfg, 'stages',scoring.label); %select everything

    %assure the stages are a cellstr structure, even with one element
    if ~iscellstr(cfg.stages)
        if ischar(cfg.stages)
            cfg.stages = {cfg.stages};
        else
            ft_error('the stages of interest must be defines as a string or a cellstr.');
        end
    end
end


%-------check for continuous data---

% data structure provided
% check if the input data is valid for this function, the input data must be raw
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'hassampleinfo', 'yes');
if isfield(data, 'trial') && isfield(data, 'time')
    % check if data structure is likely continous data
    if ((numel(data.trial) ~= 1) || ~all(size(data.sampleinfo) == [1 2])) && ~(nargout > 1)
        ft_error('data structure does not look like continuous data and has more than one trial')
    end
end
fsample = data.fsample;


if hasScoring

    sleepStagesOfInterst = cfg.stages;
    epochs = scoring.epochs;
    excluded = scoring.excluded;

    %locate data samples
    epochLengthSamples = scoring.epochlength*fsample;
    offsetSamples = scoring.dataoffset*fsample;

    hypnEpochs = 1:numel(epochs);
    hypnEpochsBeginsSamples = (((hypnEpochs - 1) * epochLengthSamples) + 1)';
    hypnEpochsEndsSamples = (hypnEpochs * epochLengthSamples)';

    %locate epochs of interest
    if istrue(cfg.usescoringexclusion)
        epochsOfInterst = hypnEpochs(ismember(epochs',sleepStagesOfInterst) & (excluded' == 0));
    else
        epochsOfInterst = hypnEpochs(ismember(epochs',sleepStagesOfInterst));
    end



    %locate samples belonging to epochs of interest
    begins = [];
    ends = [];

    if ~isempty(epochsOfInterst)
        [consecBegins, consecEnds] = consecutiveBeginsAndEnds(epochsOfInterst,1);
        begins = hypnEpochsBeginsSamples(consecBegins);
        ends = hypnEpochsEndsSamples(consecEnds);
    else
        ft_warning('requested sleep stages are not present in the scoring, will return empty begins_epochs ends_epochs.');

    end
    cfg.contbegsample = begins+offsetSamples;
    cfg.contendsample = ends+offsetSamples;

else %if no scoring, convert custom contbegtime to contbegsample

    if isfield(cfg,'contbegtime') && isfield(cfg,'contendtime')
        cfg.contbegsample=st_times2samples(data,cfg.contbegtime);
        cfg.contendsample=st_times2samples(data,cfg.contendtime);

    end

end

%select
data = ft_redefinetrial(cfg, data);

%calculate  things to provide some feedback
data_sampleinfo=data.sampleinfo;

data_trials=size(data_sampleinfo,1);
data_minutes=sum(diff(data_sampleinfo,[],2))/data.fsample/60;

fprintf('selected %i trials with minimum trial length of %i s and combined duration of %.1f min\n',data_trials,cfg.minlength,data_minutes)

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)




