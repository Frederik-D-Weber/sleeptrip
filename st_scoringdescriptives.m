function [res] = st_scoringdescriptives(cfg, scoring)
% 
% ST_SCORINGDESCRIPTIVES determine sleep scoring values for sleep table and stores 
% in result structure
%
% Use as
%   [res] = st_scoringdescriptives(cfg, scoring)
%
% Configuration parameter can be empty, e.g. cfg = []
%
% Optional configuration parameters are
%   cfg.sleeponsetdef  = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                        'NR' or 'N2R' or 'XR' or 'AASM' or 'X2R' or 
%                        'N2' or 'N3' or 'SWS' or 'S4' or 'R',
%                        see ST_SLEEPONSET for details (default = 'N1_XR')
%
% See also ST_READ_SCORING, ST_SLEEPONSET

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

tic
memtic
st = dbstack;
functionname = st.name;
fprintf([functionname ' function started\n']);


% set the defaults
cfg.sleeponsetdef  = ft_getopt(cfg, 'sleeponsetdef', 'N1_XR');


% assumes there is at least
% one valid epoch of N2 or N3 or S4 or R

hasLightsOff = false;
lightsOffMoment = 0;
if isfield(scoring, 'lightsoff')
    hasLightsOff = true;
    lightsOffMoment = scoring.lightsoff;
else
    ft_warning('The lights off moment was not provided in the scoring structure.\n The beginning of the scoring is thus assumed as lights off.');
end

sleepOnsetTime = NaN;

N1OnsetTime = NaN;
N2OnsetTime = NaN;

SWSonsetTime = NaN;
S4onsetTime = NaN;
REMonsetTime = NaN;

totalSleepTime = NaN;

N1Time = NaN;
N2Time = NaN;
N3Time = NaN;
S4Time = NaN;
REMtime = NaN;
WakeTime = NaN;
MovementTime = NaN;
ArtifactTime = NaN;
UnknownTime = NaN;

SWStime = NaN;
NRtime = NaN;

N1TimePreOnset = NaN;
N2TimePreOnset = NaN;
N3TimePreOnset = NaN;
S4TimePreOnset = NaN;
REMtimePreOnset = NaN;
WakeTimePreOnset = NaN;
MovementTimePreOnset = NaN;
ArtifactTimePreOnset = NaN;
UnknownTimePreOnset = NaN;
SWStimePreOnset = NaN;
NRtimePreOnset = NaN;

N1Time_without_excluded = NaN;
N2Time_without_excluded = NaN;
N3Time_without_excluded = NaN;
S4Time_without_excluded = NaN;
REMtime_without_excluded = NaN;
WakeTime_without_excluded = NaN;
MovementTime_without_excluded = NaN;
ArtifactTime_without_excluded = NaN;
UnknownTime_without_excluded = NaN;
SWStime_without_excluded = NaN;
NRtime_without_excluded = NaN;

LongestWakeTimePeriodAfterSleepOnset = NaN;
LongestWakeTimePeriodAfterSleepOnset_after_so = NaN;



fprintf([functionname ' function initialized\n']);


dummySampleRate = 100;
lightsOffSample = round(lightsOffMoment*dummySampleRate);
epochLengthSamples = scoring.epochlength * dummySampleRate;

epochs = cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0);
hypnStages = [cellfun(@sleepStage2str,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt2,epochs,'UniformOutput',0)];


hypnEpochs = 1:numel(epochs);
hypnEpochsBeginsSamples = (((hypnEpochs - 1) * epochLengthSamples) + 1)';
hypnEpochsEndsSamples   = (hypnEpochs * epochLengthSamples)';

[onsetCandidate preOffsetCandidate onsetepoch] = st_sleeponset(cfg,scoring);

sleepOnsetTime = (onsetCandidate-1)*scoring.epochlength - lightsOffMoment;

N1ind = find(strcmp(hypnStages(:,1),'N1'));
N2ind = find(strcmp(hypnStages(:,1),'N2'));
SWSind = find(strcmp(hypnStages(:,2),'SWS'));
S4ind = find(strcmp(hypnStages(:,1),'S4'));
REMind = find(strcmp(hypnStages(:,1),'R'));

if ~isempty(N1ind)
    N1OnsetTime = (min(N1ind(N1ind >= onsetCandidate)) - onsetCandidate)*scoring.epochlength;
    if isempty(N1OnsetTime)
        N1OnsetTime = NaN;
    end
end
if ~isempty(N2ind)
    N2OnsetTime = (min(N2ind(N2ind >= onsetCandidate)) - onsetCandidate)*scoring.epochlength;
    if isempty(N2OnsetTime)
        N2OnsetTime = NaN;
    end
end
if ~isempty(SWSind)
    SWSonsetTime = (min(SWSind(SWSind >= onsetCandidate)) - onsetCandidate)*scoring.epochlength;
    if isempty(SWSind)
        SWSind = NaN;
    end
end
if ~isempty(S4ind)
    S4onsetTime = (min(S4ind(S4ind >= onsetCandidate)) - onsetCandidate)*scoring.epochlength;
    if isempty(S4onsetTime)
        S4onsetTime = NaN;
    end
end
if ~isempty(REMind)
    REMonsetTime = (min(REMind(REMind >= onsetCandidate)) - onsetCandidate)*scoring.epochlength;
    if isempty(REMonsetTime)
        REMonsetTime = NaN;
    end
end

preOnsetCandidate = onsetCandidate;
if preOnsetCandidate > 1
    preOnsetCandidate = preOnsetCandidate-1;
end
%hypnTST = epochs(onsetCandidate:preOffsetCandidate);
hypnTST_excluded = scoring.excluded(onsetCandidate:preOffsetCandidate);


hypnStagesTST = hypnStages(onsetCandidate:preOffsetCandidate,:);
hypnStagesPreSleepOnset = hypnStages(1:preOnsetCandidate,:);
hypnEpochsTST = hypnEpochs(onsetCandidate:preOffsetCandidate);
hypnEpochsBeginsSamplesTST = hypnEpochsBeginsSamples(onsetCandidate:preOffsetCandidate,:);
hypnEpochsEndsSamplesTST = hypnEpochsEndsSamples(onsetCandidate:preOffsetCandidate,:);

totalSleepTime = (length(onsetCandidate:preOffsetCandidate))*scoring.epochlength;

N1Time = length(find(strcmp(hypnStagesTST(:,1),'N1')))*scoring.epochlength;
N2Time = length(find(strcmp(hypnStagesTST(:,1),'N2')))*scoring.epochlength;
N3Time = length(find(strcmp(hypnStagesTST(:,1),'N3')))*scoring.epochlength;
S4Time = length(find(strcmp(hypnStagesTST(:,1),'S4')))*scoring.epochlength;
REMtime = length(find(strcmp(hypnStagesTST(:,1),'R')))*scoring.epochlength;
WakeTime = length(find(strcmp(hypnStagesTST(:,1),'W')))*scoring.epochlength;
MovementTime = length(find(strcmp(hypnStagesTST(:,1),'MT')))*scoring.epochlength;
ArtifactTime = length(find(strcmp(hypnStagesTST(:,1),'A')))*scoring.epochlength;
UnknownTime = length(find(strcmp(hypnStagesTST(:,1),'?')))*scoring.epochlength;


SWStime = N3Time + S4Time;
NRtime = SWStime + N2Time;


N1TimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'N1')))*scoring.epochlength;
N2TimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'N2')))*scoring.epochlength;
N3TimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'N3')))*scoring.epochlength;
S4TimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S4')))*scoring.epochlength;
REMtimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'R')))*scoring.epochlength;
WakeTimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'W')))*scoring.epochlength;
MovementTimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'MT')))*scoring.epochlength;
ArtifactTimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'A')))*scoring.epochlength;
UnknownTimePreOnset = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'?')))*scoring.epochlength;

SWStimePreOnset = N3TimePreOnset + S4TimePreOnset;
NRtimePreOnset = SWStimePreOnset + N2TimePreOnset;

N1Time_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'N1') & (hypnTST_excluded == 0)' ))*scoring.epochlength;
N2Time_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'N2') & (hypnTST_excluded == 0)' ))*scoring.epochlength;
N3Time_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'N3') & (hypnTST_excluded == 0)' ))*scoring.epochlength;
S4Time_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'S4') & (hypnTST_excluded == 0)' ))*scoring.epochlength;
REMtime_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'R') & (hypnTST_excluded == 0)' ))*scoring.epochlength;
WakeTime_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'W') & (hypnTST_excluded == 0)' ))*scoring.epochlength;
MovementTime_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'MT') & (hypnTST_excluded == 0)' ))*scoring.epochlength;
ArtifactTime_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'A') & (hypnTST_excluded == 0)' ))*scoring.epochlength;
UnknownTime_without_excluded = length(find(strcmp(hypnStagesTST(:,1),'?') & (hypnTST_excluded == 0)' ))*scoring.epochlength;

SWStime_without_excluded = N3Time_without_excluded + S4Time_without_excluded;
NRtime_without_excluded = SWStime_without_excluded + N2Time_without_excluded;

hypnStagesTST_wake_index = find(strcmp(hypnStagesTST(:,1),'W'));
if ~isempty(hypnStagesTST_wake_index)
    [hypnStagesTST_wake_index_begins, hypnStagesTST_wake_index_ends] = consecutiveBeginsAndEnds(hypnStagesTST_wake_index,1);
    hypnStagesTST_wake_consec_epoch_duration = hypnStagesTST_wake_index_ends - hypnStagesTST_wake_index_begins + 1;
    LongestWakeTimePeriodAfterSleepOnset = max(hypnStagesTST_wake_consec_epoch_duration)*scoring.epochlength;
    LongestWakeTimePeriodAfterSleepOnset_after_so = hypnStagesTST_wake_index_begins(find(hypnStagesTST_wake_consec_epoch_duration == max(hypnStagesTST_wake_consec_epoch_duration),1,'first'))*scoring.epochlength;
else
    LongestWakeTimePeriodAfterSleepOnset = 0;
    LongestWakeTimePeriodAfterSleepOnset_after_so = NaN;
end

res = [];
st = dbstack;
res.ori = st.name;
res.type = 'descriptive';
res.cfg = cfg;
res.table = table(...
    {scoring.ori.scoringfile}, ...
    {scoring.standard}, ...
    scoring.epochlength, ...
    totalSleepTime/60, ...
    sleepOnsetTime/60, ...
    N1OnsetTime/60, ...
    N2OnsetTime/60, ...
    SWSonsetTime/60, ...
    S4onsetTime/60, ...
    REMonsetTime/60, ...
    N1Time/60, ...
    N2Time/60, ...
    N3Time/60, ...
    S4Time/60, ...
    REMtime/60, ...
    WakeTime/60, ...
    MovementTime/60, ...
    ArtifactTime/60, ...
    UnknownTime/60, ...
    SWStime/60, ...
    NRtime/60, ...
    100*N1Time/totalSleepTime, ...
    100*N2Time/totalSleepTime, ...
    100*N3Time/totalSleepTime, ...
    100*S4Time/totalSleepTime, ...
    100*REMtime/totalSleepTime, ...
    100*WakeTime/totalSleepTime, ...
    100*MovementTime/totalSleepTime, ...
    100*ArtifactTime/totalSleepTime, ...
    100*UnknownTime/totalSleepTime, ...
    100*SWStime/totalSleepTime, ...
    100*NRtime/totalSleepTime, ...
    N1Time_without_excluded/60, ...
    N2Time_without_excluded/60, ...
    N3Time_without_excluded/60, ...
    S4Time_without_excluded/60, ...
    REMtime_without_excluded/60, ...
    WakeTime_without_excluded/60, ...
    MovementTime_without_excluded/60, ...
    ArtifactTime_without_excluded/60, ...
    UnknownTime_without_excluded/60, ...
    SWStime_without_excluded/60, ...
    NRtime_without_excluded/60, ...
    100*N1Time_without_excluded/totalSleepTime, ...
    100*N2Time_without_excluded/totalSleepTime, ...
    100*N3Time_without_excluded/totalSleepTime, ...
    100*S4Time_without_excluded/totalSleepTime, ...
    100*REMtime_without_excluded/totalSleepTime, ...
    100*WakeTime_without_excluded/totalSleepTime, ...
    100*MovementTime_without_excluded/totalSleepTime, ...
    100*ArtifactTime_without_excluded/totalSleepTime, ...
    100*UnknownTime_without_excluded/totalSleepTime, ...
    100*SWStime_without_excluded/totalSleepTime, ...
    100*NRtime_without_excluded/totalSleepTime, ...
    N1TimePreOnset/60, ...
    N2TimePreOnset/60, ...
    N3TimePreOnset/60, ...
    S4TimePreOnset/60, ...
    REMtimePreOnset/60, ...
    WakeTimePreOnset/60, ...
    MovementTimePreOnset/60, ...
    ArtifactTimePreOnset/60, ...
    UnknownTimePreOnset/60, ...
    SWStimePreOnset/60, ...
    NRtimePreOnset/60, ...
    LongestWakeTimePeriodAfterSleepOnset/60, ...
    LongestWakeTimePeriodAfterSleepOnset_after_so/60, ...
    'VariableNames',{'scoringfile','standard','epoch_length_seconds','Total_sleep_time_min','Sleep_Onset_min','N1_onset_min','N2_onset_min','SWS_onset_min','S4_onset_min','R_onset_min'...
    ,'N1_min','N2_min','N3_min','S4_min','R_min','Wake_after_sleep_onset_min','Movement_Time_min','Artifact_min','Unknown_min','SWS_min','NR_without_N1_min'...
    ,'N1_percent','N2_percent','N3_percent','S4_percent','R_percent','Wake_after_sleep_onset_percent','Movement_Time_percent','Artifact_percent','Unknown_percent','SWS_percent','NR_without_N1_percent'...
    ,'N1_without_excluded_min','N2_without_excluded_min','N3_without_excluded_min','S4_without_excluded_min','R_without_excluded_min','Wake_after_sleep_onset_without_excluded_min','Movement_Time_without_excluded_min', 'Artifact_without_excluded_min', 'Unknown_without_excluded_min', 'SWS_without_excluded_min','NR_without_N1_without_excluded_min'...
    ,'N1_without_excluded_percent','N2_without_excluded_percent','N3_without_excluded_percent','S4_without_excluded_percent','R_without_excluded_percent','Wake_after_sleep_onset_without_excluded_percent','Movement_Time_without_excluded_percent', 'Artifact_without_excluded_percent', 'Unknown_without_excluded_percent','SWS_without_excluded_percent','NR_without_N1_without_excluded_percent'...
    ,'N1_before_sleep_onset_min','N2_before_sleep_onset_min','N3_before_sleep_onset_min','S4_before_sleep_onset_min','R_before_sleep_onset_min','Wake_before_sleep_onset_min','Movement_Time_before_sleep_onset_min', 'Artifact_before_sleep_onset_min', 'Unknown_before_sleep_onset_min','SWS_before_sleep_onset_min','NR_before_sleep_onset_without_N1_min'...
    ,'longest_WASO_period_min','longest_WASO_period_after_latency_after_so_min'}...
    );


fprintf([functionname ' function finished\n']);

toc
memtoc
end
