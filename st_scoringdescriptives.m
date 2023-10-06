function [res] = st_scoringdescriptives(cfg, scoring)
%
% ST_SCORINGDESCRIPTIVES determine sleep scoring values for sleep table and stores
% in result structure
% Note that many other sleep variables can be derived from these
% descriptives (e.g. sleep efficiency).
% If scoring.lightsoff and scoring.lightson fields are present (in seconds since recording onset)
% they will be used for sleep onset (latency)
% the columns N1_delay_min, N2_delay_min, SWS_delay_min, S4_delay_min,
% R_delay_min are considered time after sleep onset, i.e. if sleep onset is
% not occuring, also those onset delays are not considered.
% Percentages (perc.) are given in respect to the sleep period (i.e. duration after
% sleep onset to sleep end) consisting of wake after sleep onset (WASO) and
% total sleep duration in the sleep period and with respect to the actual sleeping time.
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
%   cfg.allowsleeponsetbeforesleepopon = srting, if possible, allow sleep onset before sleep
%                        opportunity (or lights off moment if former is not present)
%                        either 'yes' or 'no' see ST_SLEEPONSET for details (default = 'no')
%   cfg.allowsleepaftersleepopoff = srting, if possible, allow sleep (offset, i.e. end of sleep) after sleep
%                        opportunity (or lights on moment if former is not present)
%                        either 'yes' or 'no' see ST_SLEEPONSET for details (default = 'yes')
%   cfg.allowsleepopoonbeforescoring = srting, if possible, allow sleep opportunity onset to be before the scoring starts, i.e. to be < 0(default = 'yes')
%   cfg.allowsleepopoffafterscoring = srting, if possible, allow sleep opportunity offset to be after the scoring ends, i.e. to be > scoring duration(default = 'yes')
%   cfg.fixindicatorstoepochlength = string, either 'no', 'start', 'end', and 'snap' to put the
%                       indicators like lightsoff, lightson, sleepopon and sleepopoff in the scoring
%                       to the start or end of their epoch or snap them to either the start or end on which whatever is closer (default = 'no')
%   cfg.cycle          = string or integer value from 1 to Inf.
%                        If set, then the descriptives is derived from
%                        only this cycle.
%                        possible values are:
%                          'presleep' = prior to first cycle
%                          0 = presleep cycle
%                          1 = first cycle in sleep
%                          2 = second cycle in slepp
%                          3 = third cycle in sleep
%                          [1 3] = first and third cycle in sleep
%                          ...
%                          'last' = last cycle (even uncomplete) in sleep%
%                          'last' = last complete cycle in sleep
%                          'postsleep' = post last cycle (after last)
%                          'all' = all sleep cycles (presleep,1,....,last,postsleep)
%                        See ST_SLEEPCYCLES for futher
%                        parameters that are passed on to it.
%
%
% See also ST_READ_SCORING, ST_SLEEPONSET, ST_SLEEPCYCLES

% Copyright (C) 2019-, Frederik D. Weber
% Thanks Roy Cox, Leonore Bovy, Shervin Bukhari for valuable suggestions
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
functionname = getfunctionname();
fprintf([functionname ' function started\n']);


% set the defaults
cfg.sleeponsetdef  = ft_getopt(cfg, 'sleeponsetdef', 'N1_XR');

cfg.allowsleepopoonbeforescoring  = ft_getopt(cfg, 'allowsleepopoonbeforescoring', 'yes');
cfg.allowsleepopoffafterscoring  = ft_getopt(cfg, 'allowsleepopoffafterscoring', 'yes');

cfg.fixindicatorstoepochlength  = ft_getopt(cfg, 'fixindicatorstoepochlength', 'no');

if ~strcmp(cfg.fixindicatorstoepochlength,'no')
    ft_warning('trying to fix the indicators like lightsoff, lightson, sleepopon and sleepopoff in the scoring to whole epochs with ''%s''.',cfg.fixindicatorstoepochlength)
    if isfield(scoring, 'lightsoff')
        switch cfg.fixindicatorstoepochlength
            case 'snap'
                scoring.lightsoff = round(scoring.lightsoff/scoring.epochlength)*scoring.epochlength;
            case 'start'
                scoring.lightsoff = floor(scoring.lightsoff/scoring.epochlength)*scoring.epochlength;
            case 'begin'
                scoring.lightsoff = ceil(scoring.lightsoff/scoring.epochlength)*scoring.epochlength;
        end
    end
    if isfield(scoring, 'lightson')
        switch cfg.fixindicatorstoepochlength
            case 'snap'
                scoring.lightson = round(scoring.lightson/scoring.epochlength)*scoring.epochlength;
            case 'start'
                scoring.lightson = floor(scoring.lightson/scoring.epochlength)*scoring.epochlength;
            case 'begin'
                scoring.lightson = ceil(scoring.lightson/scoring.epochlength)*scoring.epochlength;
        end
    end
    if isfield(scoring, 'sleepopon')
        switch cfg.fixindicatorstoepochlength
            case 'snap'
                scoring.sleepopon = round(scoring.sleepopon/scoring.epochlength)*scoring.epochlength;
            case 'start'
                scoring.sleepopon = floor(scoring.sleepopon/scoring.epochlength)*scoring.epochlength;
            case 'begin'
                scoring.sleepopon = ceil(scoring.sleepopon/scoring.epochlength)*scoring.epochlength;
        end
    end
    if isfield(scoring, 'sleepopoff')
        switch cfg.fixindicatorstoepochlength
            case 'snap'
                scoring.sleepopoff = round(scoring.sleepopoff/scoring.epochlength)*scoring.epochlength;
            case 'start'
                scoring.sleepopoff = floor(scoring.sleepopoff/scoring.epochlength)*scoring.epochlength;
            case 'begin'
                scoring.sleepopoff = ceil(scoring.sleepopoff/scoring.epochlength)*scoring.epochlength;
        end
    end
end


% assumes there is at least
% one valid epoch of N2 or N3 or S4 or R

hasLightsOff = false;
lightsOffMoment = NaN;
if isfield(scoring, 'lightsoff')
    if ~isnan(scoring.lightsoff)
        hasLightsOff = true;
        lightsOffMoment = scoring.lightsoff;
    else
        ft_warning('The lights off moment was NaN in the scoring structure.\n The beginning of the scoring is thus assumed as lights off.');
    end
else
    ft_warning('The lights off moment was not provided in the scoring structure.\n The beginning of the scoring is thus assumed as lights off.');
end

hasSleepOpportunityOn = false;
sleepOpportunityOnMoment = NaN;
if isfield(scoring, 'sleepopon')
    if ~isnan(scoring.sleepopon)
        hasSleepOpportunityOn = true;
        sleepOpportunityOnMoment = scoring.sleepopon;
    else
        ft_warning('The sleep opportunity onset moment was NaN in the scoring structure.\n The beginning of the scoring is thus assumed as sleep opportunity onset, but sleep onset latency will be NaN.');
    end
else
    if hasLightsOff && ~isnan(lightsOffMoment)
        sleepOpportunityOnMoment = lightsOffMoment;
        ft_warning('The sleep opportunity onset moment was not provided in the scoring structure.\n The lights off moment is used instead');
    else
        ft_warning('The sleep opportunity onset moment was not provided in the scoring structure.\n The beginning of the scoring is thus assumed as sleep opportunity onset, but sleep onset latency will be NaN.');
    end
end

scoring_cycles = {scoring};
cycle_complete = NaN;
cycle_labels   = {'non-separated'};

scoring_duration_min = numel(scoring.epochs)*scoring.epochlength/60;

data_offset_min = NaN;
if isfield(scoring,'dataoffset')
    data_offset_min = scoring.dataoffset/60;
end




if isfield(cfg,'cycle')
    res_cycle = st_sleepcycles(cfg,scoring);
    cycle_complete = res_cycle.table.complete;
    cycle_labels = arrayfun(@(x) cellstr(num2str(x)),res_cycle.table.cycle);
    cfg2 = [];
    cfg2.start = res_cycle.table.startepoch;
    cfg2.end = res_cycle.table.endepoch;

    if isempty(cfg2.start)
        ft_warning('No sleep cycles detected in this scoring')
    else
        if (cfg2.start(1) > 1)
            cfg2.end = [(cfg2.start(1) - 1); cfg2.end];
            cfg2.start = [1 ; cfg2.start];
            cycle_complete = [false  ; cycle_complete];
            cycle_labels = ['presleep' ; cycle_labels];
        end
        if (cfg2.end(end) < numel(scoring.epochs))
            cfg2.start = [cfg2.start ; (cfg2.end(end) + 1)];
            cfg2.end = [cfg2.end ; numel(scoring.epochs)];
            cycle_complete = [cycle_complete ; false];
            cycle_labels = [cycle_labels ; 'postsleep'];
        end
    end
    scoring_cycles = st_cutscoring(cfg2,scoring);

    scoring_duration_min = (cfg2.end - cfg2.start + 1)*scoring.epochlength/60;

    indcycle = [];
    if isnumeric(cfg.cycle)
        if numel(cfg.cycle) > 1
            if any(strcmp(cycle_labels,'presleep'))
                indcycle = cfg.cycle+1;
            else
                indcycle = cfg.cycle;
            end
        else
            if cfg.cycle == 0
                indcycle = strcmp(cycle_labels,'presleep');
            else
                indcycle = ismember(cycle_labels,arrayfun(@num2str,cfg.cycle,'UniformOutput',false));
            end
        end
    else
        switch cfg.cycle
            case 'last'
                indcycle = logical(zeros(numel(cycle_complete),1));
                indcycle(find(~cellfun(@isempty,(cellfun(@str2num,cycle_labels,'UniformOutput',false))),1,'last')) = true;
            case 'lastcomplete'
                indcycle = logical(zeros(numel(cycle_complete),1));
                indcycle(find(cycle_complete,1,'last')) = true;
            case 'presleep'
                indcycle = strcmp(cycle_labels,'presleep');
            case 'postsleep'
                indcycle = strcmp(cycle_labels,'postsleep');
            case 'complete'
                indcycle = cycle_complete;
            case 'all'
                indcycle = 1:numel(scoring_cycles);
            otherwise
                ft_error('cfg.cycle is unknown, please redefine')
        end
    end


    data_offset_min_cycle = nan(1,numel(cfg2.end));
    for iSc = 1:numel(cfg2.end)
        if ~isempty(scoring_duration_min)
            data_offset_min_cycle(iSc) = data_offset_min + sum(scoring_duration_min(1:(iSc-1)));
        end
    end

    scoring_cycles = scoring_cycles(indcycle);
    cycle_complete = cycle_complete(indcycle);
    cycle_labels = cycle_labels(indcycle);
    scoring_duration_min = scoring_duration_min(indcycle);
    data_offset_min = data_offset_min_cycle(indcycle);
end


res = [];
res.ori = getfunctionname();
res.type = 'descriptive';
res.cfg = cfg;

for iScoringCycle = 1:numel(scoring_cycles)
    scoring = scoring_cycles{iScoringCycle};
    res_cycle_sccy = st_sleepcycles(cfg,scoring);
    if isempty(res_cycle_sccy.table.startepoch)
        ft_warning('No sleep cycles detected in this scoring part %d',iScoringCycle)
    end
    sleepOnsetTime = NaN;

    N1OnsetTime = NaN;
    N2OnsetTime = NaN;

    SWSonsetTime = NaN;
    S4onsetTime = NaN;
    REMonsetTime = NaN;

    totalSleepPeriod = NaN;

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


    N1Time_all = NaN;
    N2Time_all = NaN;
    N3Time_all = NaN;
    S4Time_all = NaN;
    REMtime_all = NaN;
    WakeTime_all = NaN;
    MovementTime_all = NaN;
    ArtifactTime_all = NaN;
    UnknownTime_all = NaN;
    SWStime_all = NaN;
    NRtime_all = NaN;

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

    nArousals = NaN;
    nArousals_sleep_period = NaN;

    stagesCombo = {{'N1','S1'}, {'N2','S2'}, {'N3','S3'}, {'S4'}, {'R'}, {'W'}};
    nArousalsByStage = nan(size(stagesCombo));
    nArousalsByStageNonMTnonExcl = nan(size(stagesCombo));


    fprintf([functionname ' function initialized\n']);


    dummySampleRate = 100;
    %lightsOffSample = round(lightsOffMoment*dummySampleRate);
    epochLengthSamples = scoring.epochlength * dummySampleRate;

    epochs = cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0);
    hypnStages = [cellfun(@sleepStage2str,epochs,'UniformOutput',0) ...
        cellfun(@sleepStage2str_alt,epochs,'UniformOutput',0) ...
        cellfun(@sleepStage2str_alt2,epochs,'UniformOutput',0) ...
        cellfun(@sleepStage2str_alt3,epochs,'UniformOutput',0)];


    N1Time_all = length(find(strcmp(hypnStages(:,1),'N1')))*scoring.epochlength;
    N2Time_all = length(find(strcmp(hypnStages(:,1),'N2')))*scoring.epochlength;
    N3Time_all = length(find(strcmp(hypnStages(:,1),'N3')))*scoring.epochlength;
    S4Time_all = length(find(strcmp(hypnStages(:,1),'S4')))*scoring.epochlength;
    REMtime_all = length(find(strcmp(hypnStages(:,1),'R')))*scoring.epochlength;
    WakeTime_all = length(find(strcmp(hypnStages(:,1),'W')))*scoring.epochlength;
    MovementTime_all = length(find(strcmp(hypnStages(:,1),'MT')))*scoring.epochlength;
    ArtifactTime_all = length(find(strcmp(hypnStages(:,1),'A')))*scoring.epochlength;
    UnknownTime_all = length(find(strcmp(hypnStages(:,1),'?')))*scoring.epochlength;

    SWStime_all = N3Time_all + S4Time_all;
    NRtime_all = SWStime_all + N2Time_all;

    hypnEpochs = 1:numel(epochs);
    hypnEpochsBeginsSamples = (((hypnEpochs - 1) * epochLengthSamples) + 1)';
    hypnEpochsEndsSamples   = (hypnEpochs * epochLengthSamples)';

    [onsetnumber lastscoredsleepstagenumber onsetepoch lastscoredsleepstage allowedsleeponsetbeforesleepopon allowedsleepaftersleepopoff] = st_sleeponset(cfg,scoring);

    if isempty(onsetepoch)
        onsetnumber = NaN;
    end

    hasLightsOn = false;
    lightsOnMoment = lastscoredsleepstagenumber*scoring.epochlength;
    lightsOnToLightsOff = NaN;
    if isfield(scoring, 'lightson')
        if ~isnan(scoring.lightson)
            hasLightsOn = true;
            lightsOnMoment = scoring.lightson;
            if ~istrue(cfg.allowsleepopoonbeforescoring)
                lightsOffMoment = 0;
            end
            if ~istrue(cfg.allowsleepopoffafterscoring)
                lightsOnMoment = numel(scoring.epochs)*scoring.epochlength;
            end
            lightsOnToLightsOff = (lightsOnMoment-lightsOffMoment);
        else
            ft_warning('The lights on moment was NaN in the scoring structure.\n The end of sleep is thus assumed as lights on.');
        end
    else
        ft_warning('The lights on moment was not provided in the scoring structure.\n The end of sleep is thus assumed as lights on.');
    end


    hasSleepOpportunityOff = false;
    sleepOpportunityOffMoment = lastscoredsleepstagenumber*scoring.epochlength;
    sleepOpportunityOnToSleepOpportunityOff = NaN;
    if isfield(scoring, 'sleepopoff')
        if ~isnan(scoring.sleepopoff)
            hasSleepOpportunityOff = true;
            sleepOpportunityOffMoment = scoring.sleepopoff;
            if ~istrue(cfg.allowsleepopoonbeforescoring)
                sleepOpportunityOnMoment = 0;
            end
            if ~istrue(cfg.allowsleepopoffafterscoring)
                sleepOpportunityOffMoment = numel(scoring.epochs)*scoring.epochlength;
            end
            sleepOpportunityOnToSleepOpportunityOff = (sleepOpportunityOffMoment-sleepOpportunityOnMoment);
        else
            ft_warning('The sleep opportunity off moment was NaN in the scoring structure.\n The end of sleep is thus assumed as sleep opportunity off.');
        end
    else
        if hasLightsOff && ~isnan(lightsOnMoment)
            sleepOpportunityOffMoment = lightsOnMoment;
            ft_warning('The sleep opportunity off moment was not provided in the scoring structure.\n The lights on moment is used instead');
        else
            ft_warning('The sleep opportunity off moment was not provided in the scoring structure.\n  The end of sleep is thus assumed as sleep opportunity off.');
        end
    end


    sleepOnsetTime = (onsetnumber-1)*scoring.epochlength - sleepOpportunityOnMoment;
    sleepOnsetTimepointafterscoringstart = (onsetnumber-1)*scoring.epochlength;
    sleepOffsetTimepointafterscoringstart = lastscoredsleepstagenumber*scoring.epochlength;


    if allowedsleeponsetbeforesleepopon
        ft_warning('sleep onset was allowed to be before sleep opportunity moment and thus sleep onset latency and sleep opportunity duration will be NaN.')
        sleepOnsetTime = NaN;
        sleepOpportunityOnToSleepOpportunityOff = NaN;
    end

    N1ind = find(strcmp(hypnStages(:,1),'N1'));
    N2ind = find(strcmp(hypnStages(:,1),'N2'));
    SWSind = find(strcmp(hypnStages(:,2),'SWS'));
    S4ind = find(strcmp(hypnStages(:,1),'S4'));
    REMind = find(strcmp(hypnStages(:,1),'R'));

    if ~isempty(N1ind)
        N1OnsetTime = (min(N1ind(N1ind >= onsetnumber)) - onsetnumber)*scoring.epochlength;
        if isempty(N1OnsetTime)
            N1OnsetTime = NaN;
        end
    end
    if ~isempty(N2ind)
        N2OnsetTime = (min(N2ind(N2ind >= onsetnumber)) - onsetnumber)*scoring.epochlength;
        if isempty(N2OnsetTime)
            N2OnsetTime = NaN;
        end
    end
    if ~isempty(SWSind)
        SWSonsetTime = (min(SWSind(SWSind >= onsetnumber)) - onsetnumber)*scoring.epochlength;
        if isempty(SWSonsetTime)
            SWSonsetTime = NaN;
        end
    end
    if ~isempty(S4ind)
        S4onsetTime = (min(S4ind(S4ind >= onsetnumber)) - onsetnumber)*scoring.epochlength;
        if isempty(S4onsetTime)
            S4onsetTime = NaN;
        end
    end
    if ~isempty(REMind)
        REMonsetTime = (min(REMind(REMind >= onsetnumber)) - onsetnumber)*scoring.epochlength;
        if isempty(REMonsetTime)
            REMonsetTime = NaN;
        end
    end

    preOnsetCandidate = onsetnumber;
    if preOnsetCandidate > 1
        preOnsetCandidate = preOnsetCandidate-1;
    end
    if ~isnan(onsetnumber)
        %hypnTST = epochs(onsetCandidate:preOffsetCandidate);
        hypnTST_excluded = scoring.excluded(onsetnumber:lastscoredsleepstagenumber);

        hypnStagesTST = hypnStages(onsetnumber:lastscoredsleepstagenumber,:);
        hypnStagesPreSleepOnset = hypnStages(1:preOnsetCandidate,:);
        hypnEpochsTST = hypnEpochs(onsetnumber:lastscoredsleepstagenumber);
        hypnEpochsBeginsSamplesTST = hypnEpochsBeginsSamples(onsetnumber:lastscoredsleepstagenumber,:);
        hypnEpochsEndsSamplesTST = hypnEpochsEndsSamples(onsetnumber:lastscoredsleepstagenumber,:);

        totalSleepPeriod = (length(onsetnumber:lastscoredsleepstagenumber))*scoring.epochlength;
    else

        %hypnTST = epochs(onsetCandidate:preOffsetCandidate);
        hypnTST_excluded = scoring.excluded(1:0);
        hypnStagesTST = hypnStages(1:0,:);
        hypnStagesPreSleepOnset = hypnStages(1:numel(scoring.epochs),:);
        hypnEpochsTST = hypnEpochs(1:0);
        hypnEpochsBeginsSamplesTST = hypnEpochsBeginsSamples(1:0,:);
        hypnEpochsEndsSamplesTST = hypnEpochsEndsSamples(1:0,:);

        totalSleepPeriod = 0;
    end


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


    if isfield(scoring, 'arousals')
        nArousals = size(scoring.arousals,1);

        for iStages = 1:numel(stagesCombo)


            if ~isnan(onsetnumber)
                cfg_cut = [];
                cfg_cut.start = onsetnumber;
                cfg_cut.end = lastscoredsleepstagenumber;
                scorings_sleepperiod = st_cutscoring(cfg_cut,scoring);
                scorings_sleepperiod = scorings_sleepperiod{1};
            else
                scorings_sleepperiod = scoring;
                scorings_sleepperiod.epochs = scorings_sleepperiod.epochs(1:0);
                scorings_sleepperiod.excluded = scorings_sleepperiod.excluded(1:0);
            end

            nArousals_sleep_period = size(scorings_sleepperiod.arousals,1);

            cfg_sc = cfg;
            cfg_sc.stages = stagesCombo{iStages};
            cfg_sc.considerexcluded = 'no';
            [begins_epoch ends_epoch] = st_select_scoring(cfg_sc,scorings_sleepperiod);
            if ~isempty(begins_epoch)
                cfg_sc.start = begins_epoch;
                cfg_sc.end = ends_epoch;
                scorings_cut = st_cutscoring(cfg_sc,scorings_sleepperiod);

                scorings_cut_appended = st_append_scoring(scorings_cut{:});
                nArousalsByStage(iStages) = size(scorings_cut_appended.arousals,1);
            end

            cfg_sc = cfg;
            cfg_sc.stages = stagesCombo{iStages};
            cfg_sc.considerexcluded = 'yes';
            [begins_epoch ends_epoch] = st_select_scoring(cfg_sc,scorings_sleepperiod);
            if ~isempty(begins_epoch)
                cfg_sc.start = begins_epoch;
                cfg_sc.end = ends_epoch;
                scorings_cut = st_cutscoring(cfg_sc,scorings_sleepperiod);

                scorings_cut_appended = st_append_scoring(scorings_cut{:});
                nArousalsByStageNonMTnonExcl(iStages) = size(scorings_cut_appended.arousals,1);
            end
        end
    end

    hypnEpochs = 1:numel(scoring.epochs);
    hypnEpochsTST = 1:numel(hypnStagesTST(:,1));


    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'N1'})),1);
    N1_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'N2'})),1);
    N2_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'N3'})),1);
    N3_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'S4'})),1);
    S4_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'R'})),1);
    REM_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'W'})),1);
    Wake_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'MT'})),1);
    MovementTime_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'A'})),1);
    Artifact_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'?'})),1);
    Unknown_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'N3','S4'})),1);
    SWS_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'N2','N3','S4'})),1);
    NR_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'N1','N2','N3','S4'})),1);
    NRwithN1_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,4),{'S'})),1);
    S_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,4),{'W'})),1);
    nonS_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochsTST(ismember(hypnStagesTST(:,1),{'MT','W'})),1);
    WakeMT_bouts = numel(consecBegins);

    idx_unknown_TST = hypnEpochsTST(ismember(hypnStagesTST(:,1),{'?'}));
    [temp1, temp2, arbitrary_numbers] = unique(hypnStagesTST(:,1));
    transtitions = logical(abs([0; diff(arbitrary_numbers)]));
    transtitions(idx_unknown_TST) = 0;

    NtransitionsTST = sum(transtitions);

    %------ sleep fragmentation-----------
    nCycles = size(res_cycle_sccy.table,1);

    interruptTypes = {{'MT','W'},{'MT','W','N1','S1'},{'MT','W','N1','S1','N2','S2','N3','S3','S4'}};
    episodeTypes_start = {'NRstartepoch', 'Rstartepoch', 'NRstarts_ascending'};
    episodeType_ends = {'NRendepoch', 'Rendepoch', 'NRends_ascending'};
    nInterruptTypes = numel(interruptTypes);
    nEpisodeTypes = numel(episodeTypes_start);

    interrupt_duration_min_by_IntType_and_EpType = cell(nInterruptTypes,nEpisodeTypes,2);
    timewindow_min_by_IntType_and_EpType = nan(nInterruptTypes,nEpisodeTypes,2);
    Interruptions_in_bouts_by_IntType_and_EpType = nan(nInterruptTypes,nEpisodeTypes,2);
    mean_Duration_min_by_IntType_and_EpType = nan(nInterruptTypes,nEpisodeTypes,2);
    overall_Duration_min_by_IntType_and_EpType = nan(nInterruptTypes,nEpisodeTypes,2);

    %for each definition of an interruption (=constellation of stages)
    for iInterruptType = 1:nInterruptTypes
        interruptType = interruptTypes{iInterruptType};


        interrupt_duration_min_byEpType_and_Cycle = cell(nEpisodeTypes,nCycles);
        overall_Duration_min_byEpType_and_Cycle = nan(nEpisodeTypes,nCycles);
        Interruptions_in_bouts_byEpType_and_Cycle = nan(nEpisodeTypes,nCycles);
        mean_Duration_min_byEpType_and_Cycle = nan(nEpisodeTypes,nCycles);
        timewindow_seconds_byEpType_and_Cycle = nan(nEpisodeTypes,nCycles);

        %for each definition of an episode (= NREM, REM, ascendingNREM)
        for iEpisodeType = 1:numel(episodeTypes_start)

            %for two variants: excluding/including manual arousals
            for iArousalNoYes = 1:2

                %skip combination of NREM episode plus NREM-interruptions
                if strcmp(episodeTypes_start{iEpisodeType},'NRstartepoch') && (iInterruptType == 3)
                    continue;
                end

                %skip combination of REM episode plus N1-interruptions (not sure why)
                if strcmp(episodeTypes_start{iEpisodeType},'Rstartepoch') && (iInterruptType == 2)
                    continue;
                end

                %for each cycle, get start and end epochs for current episode
                rangestart = res_cycle_sccy.table.(episodeTypes_start{iEpisodeType});
                rangestop  = res_cycle_sccy.table.(episodeType_ends{iEpisodeType});
                for iCycle = 1:nCycles

                    scoringpart_range = rangestart(iCycle):rangestop(iCycle);

                    %scoringpart_range may be NaN. only use numel if it isn't
                    if ~isnan(scoringpart_range)
                        timewindow_seconds_byEpType_and_Cycle(iEpisodeType,iCycle) = numel(scoringpart_range)*scoring.epochlength;
                    else
                        timewindow_seconds_byEpType_and_Cycle(iEpisodeType,iCycle)=0;
                    end

                    if isnan(rangestart(iCycle)) || isnan(rangestop(iCycle)) || isempty(scoringpart_range)
                        overall_Duration_min_byEpType_and_Cycle(iEpisodeType,iCycle) = NaN;
                        Interruptions_in_bouts_byEpType_and_Cycle(iEpisodeType,iCycle) = NaN;
                        mean_Duration_min_byEpType_and_Cycle(iEpisodeType,iCycle) = NaN;
                    else

                        % if consider arousal an interruption too
                        nArousals_current_cycle = NaN;
                        if isfield(scoring, 'arousals') && (iArousalNoYes == 2)
                            cfg_sc2 = cfg;
                            cfg_sc2.considerexcluded = 'no';
                            cfg_sc2.start = scoringpart_range(1);
                            cfg_sc2.end = scoringpart_range(end);
                            scorings_cut2 = st_cutscoring(cfg_sc2,scoring);
                            scorings_cut2_appended = st_append_scoring(scorings_cut2{:});
                            nArousals_current_cycle = size(scorings_cut2_appended.arousals,1);
                        end

                        %find the subrange(s) of interruptions, and their number and total duration
                        episodeStages=hypnStages(scoringpart_range,1); %stage labels of current episode
                        [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(episodeStages,interruptType)),1);

                        consecDur = consecEnds' - consecBegins' + 1;
                        nInterrupts = numel(consecBegins);

                        if isempty(consecBegins) && isnan(nArousals_current_cycle)
                            interupptions_durations_min = [];
                            interrupt_duration_min_byEpType_and_Cycle{iEpisodeType,iCycle} = interupptions_durations_min;
                        elseif isempty(consecBegins) && ~isnan(nArousals_current_cycle) % just arousals but no other interruptions
                            if nArousals_current_cycle > 0
                                interupptions_durations_min = scorings_cut2_appended.arousals.duration/60;
                            else
                                interupptions_durations_min = [];
                            end

                        else % arousals and other interruptions


                            consecBegins_wholescoring = consecBegins + min(scoringpart_range)-1;
                            consecEnds_wholescoring = consecEnds + min(scoringpart_range)-1;

                            consecBegins_wholescoring_seconds = (consecBegins_wholescoring'-1)*scoring.epochlength;
                            consecEnds_wholescoring_seconds = consecEnds_wholescoring'*scoring.epochlength;

                            if isfield(scoring, 'arousals') && (iArousalNoYes == 2)

                                scorings_cut2_appended_temp = scorings_cut2_appended;
                                for iInterrupt = 1:numel(consecBegins_wholescoring_seconds)
                                    scorings_cut2_appended_temp.arousals(((scorings_cut2_appended_temp.arousals.start >= consecBegins_wholescoring_seconds(iInterrupt) ) & (scorings_cut2_appended_temp.arousals.stop <= consecEnds_wholescoring_seconds(iInterrupt) )),:) = [];
                                    idx_rightsideoverlap = (scorings_cut2_appended_temp.arousals.start <= consecEnds_wholescoring_seconds(iInterrupt)) & (scorings_cut2_appended_temp.arousals.start >= consecBegins_wholescoring_seconds(iInterrupt)) & (scorings_cut2_appended_temp.arousals.stop >= consecEnds_wholescoring_seconds(iInterrupt));
                                    idx_leftsideoverlap = (scorings_cut2_appended_temp.arousals.stop > consecBegins_wholescoring_seconds(iInterrupt)) & (scorings_cut2_appended_temp.arousals.stop <= consecEnds_wholescoring_seconds(iInterrupt)) & (scorings_cut2_appended_temp.arousals.start < consecBegins_wholescoring_seconds(iInterrupt));

                                    if sum(idx_rightsideoverlap) > 0
                                        consecEnds_wholescoring_seconds(iInterrupt) = max(scorings_cut2_appended_temp.arousals.stop(idx_rightsideoverlap));
                                        scorings_cut2_appended_temp.arousals(idx_rightsideoverlap,:) = [];
                                    end
                                    if sum(idx_leftsideoverlap) > 0
                                        consecBegins_wholescoring_seconds(iInterrupt) = min(scorings_cut2_appended_temp.arousals.start(idx_leftsideoverlap));
                                        scorings_cut2_appended_temp.arousals(idx_leftsideoverlap,:) = [];

                                    end
                                    scorings_cut2_appended_temp.arousals.stop(idx_leftsideoverlap) = consecBegins_wholescoring_seconds(iInterrupt);
                                end

                                consecDur_wholescoring_seconds  = consecEnds_wholescoring_seconds - consecBegins_wholescoring_seconds;

                                interupptions_durations_min = cat(1,consecDur_wholescoring_seconds/60,scorings_cut2_appended_temp.arousals.duration/60);

                            else
                                consecDur_wholescoring_seconds  = consecEnds_wholescoring_seconds - consecBegins_wholescoring_seconds;
                                interupptions_durations_min = cat(1,consecDur_wholescoring_seconds/60);
                            end
                        end

                        overall_Duration_min = (nansum(interupptions_durations_min)); %combined length of interruptions in current bout
                        Interruptions_in_bouts = numel(interupptions_durations_min(~isnan(interupptions_durations_min)));
                        if Interruptions_in_bouts >0
                            mean_Duration_min = overall_Duration_min/Interruptions_in_bouts;
                        else
                            mean_Duration_min = 0; %RC: should be probably NaN?
                            Interruptions_in_bouts = 0;
                        end
                        overall_Duration_min_byEpType_and_Cycle(iEpisodeType,iCycle) = overall_Duration_min;
                        Interruptions_in_bouts_byEpType_and_Cycle(iEpisodeType,iCycle) = Interruptions_in_bouts;
                        mean_Duration_min_byEpType_and_Cycle(iEpisodeType,iCycle) = mean_Duration_min;
                        interrupt_duration_min_byEpType_and_Cycle{iEpisodeType,iCycle} = interupptions_durations_min;
                    end
                end

                interrupt_duration_min_by_IntType_and_EpType{iInterruptType,iEpisodeType,iArousalNoYes} = cat(1,interrupt_duration_min_byEpType_and_Cycle{iEpisodeType,:});

                %aggregate across cycles, only if non-NaNs available (else nansum/nanmean changes NaN to zero)
                if ~all(isnan(Interruptions_in_bouts_byEpType_and_Cycle(iEpisodeType,:)))
                    overall_Duration_min_by_IntType_and_EpType(iInterruptType,iEpisodeType,iArousalNoYes) = nansum(interrupt_duration_min_by_IntType_and_EpType{iInterruptType,iEpisodeType,iArousalNoYes});
                    mean_Duration_min_by_IntType_and_EpType(iInterruptType,iEpisodeType,iArousalNoYes) = nanmean(interrupt_duration_min_by_IntType_and_EpType{iInterruptType,iEpisodeType,iArousalNoYes});
                    Interruptions_in_bouts_by_IntType_and_EpType(iInterruptType,iEpisodeType,iArousalNoYes) = nansum(Interruptions_in_bouts_byEpType_and_Cycle(iEpisodeType,:));
                    timewindow_min_by_IntType_and_EpType(iInterruptType,iEpisodeType,iArousalNoYes) = nansum(timewindow_seconds_byEpType_and_Cycle(iEpisodeType,:))/60;

                end
            end

        end


    end


    %---presleep------
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

    scoringfile = 'unknown';
    if isfield(scoring,'ori')
        if isfield(scoring.ori,'scoringfile')
            scoringfile = scoring.ori.scoringfile;
        end
    end


    totalSleepInSleepPeriod = (totalSleepPeriod-WakeTime);
    totalSleepInSleepPeriod_minus_unknown_or_artifact = totalSleepInSleepPeriod-ArtifactTime-UnknownTime;

    sleepEfficiency_SP_percent = 100*(totalSleepInSleepPeriod/totalSleepPeriod);
    sleepEfficiency_SOP_percent = 100*(totalSleepInSleepPeriod/(sleepOpportunityOnToSleepOpportunityOff));

    table_tmp = table(...
        {scoringfile}, ...
        {scoring.standard}, ...
        {cfg.sleeponsetdef}, ...
        scoring.epochlength, ...
        data_offset_min(iScoringCycle), ...
        scoring_duration_min(iScoringCycle), ...
        sleepOpportunityOnToSleepOpportunityOff/60, ...
        lightsOnToLightsOff/60, ...
        totalSleepPeriod/60, ...
        totalSleepInSleepPeriod/60, ...
        sleepEfficiency_SP_percent, ...
        sleepEfficiency_SOP_percent, ...
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
        (NRtime+N1Time)/60, ...
        NRtime/60, ...
        N1Time_all/60, ...
        N2Time_all/60, ...
        N3Time_all/60, ...
        S4Time_all/60, ...
        REMtime_all/60, ...
        WakeTime_all/60, ...
        MovementTime_all/60, ...
        ArtifactTime_all/60, ...
        UnknownTime_all/60, ...
        SWStime_all/60, ...
        (NRtime_all+N1Time_all)/60, ...
        NRtime_all/60, ...
        100*N1Time/totalSleepPeriod, ...
        100*N2Time/totalSleepPeriod, ...
        100*N3Time/totalSleepPeriod, ...
        100*S4Time/totalSleepPeriod, ...
        100*REMtime/totalSleepPeriod, ...
        100*WakeTime/totalSleepPeriod, ...
        100*MovementTime/totalSleepPeriod, ...
        100*ArtifactTime/totalSleepPeriod, ...
        100*UnknownTime/totalSleepPeriod, ...
        100*SWStime/totalSleepPeriod, ...
        100*(NRtime+N1Time)/totalSleepPeriod, ...
        100*NRtime/totalSleepPeriod, ...
        100*N1Time/totalSleepInSleepPeriod, ...
        100*N2Time/totalSleepInSleepPeriod, ...
        100*N3Time/totalSleepInSleepPeriod, ...
        100*S4Time/totalSleepInSleepPeriod, ...
        100*REMtime/totalSleepInSleepPeriod, ...
        100*WakeTime/totalSleepInSleepPeriod, ...
        100*MovementTime/totalSleepInSleepPeriod, ...
        100*ArtifactTime/totalSleepInSleepPeriod, ...
        100*UnknownTime/totalSleepInSleepPeriod, ...
        100*SWStime/totalSleepInSleepPeriod, ...
        100*(NRtime+N1Time)/totalSleepInSleepPeriod, ...
        100*NRtime/totalSleepInSleepPeriod, ...
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
        (NRtime_without_excluded+N1Time_without_excluded)/60, ...
        NRtime_without_excluded/60, ...
        100*N1Time_without_excluded/totalSleepPeriod, ...
        100*N2Time_without_excluded/totalSleepPeriod, ...
        100*N3Time_without_excluded/totalSleepPeriod, ...
        100*S4Time_without_excluded/totalSleepPeriod, ...
        100*REMtime_without_excluded/totalSleepPeriod, ...
        100*WakeTime_without_excluded/totalSleepPeriod, ...
        100*MovementTime_without_excluded/totalSleepPeriod, ...
        100*ArtifactTime_without_excluded/totalSleepPeriod, ...
        100*UnknownTime_without_excluded/totalSleepPeriod, ...
        100*SWStime_without_excluded/totalSleepPeriod, ...
        100*(NRtime_without_excluded+N1Time_without_excluded)/totalSleepPeriod, ...
        100*NRtime_without_excluded/totalSleepPeriod, ...
        100*N1Time_without_excluded/totalSleepInSleepPeriod, ...
        100*N2Time_without_excluded/totalSleepInSleepPeriod, ...
        100*N3Time_without_excluded/totalSleepInSleepPeriod, ...
        100*S4Time_without_excluded/totalSleepInSleepPeriod, ...
        100*REMtime_without_excluded/totalSleepInSleepPeriod, ...
        100*WakeTime_without_excluded/totalSleepInSleepPeriod, ...
        100*MovementTime_without_excluded/totalSleepInSleepPeriod, ...
        100*ArtifactTime_without_excluded/totalSleepInSleepPeriod, ...
        100*UnknownTime_without_excluded/totalSleepInSleepPeriod, ...
        100*SWStime_without_excluded/totalSleepInSleepPeriod, ...
        100*(NRtime_without_excluded+N1Time_without_excluded)/totalSleepInSleepPeriod, ...
        100*NRtime_without_excluded/totalSleepInSleepPeriod, ...
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
        allowedsleeponsetbeforesleepopon, ...
        allowedsleepaftersleepopoff, ...
        hasSleepOpportunityOn, ...
        hasSleepOpportunityOff, ...
        hasLightsOff, ...
        hasLightsOn, ...
        sleepOpportunityOnMoment/60, ...
        sleepOpportunityOffMoment/60, ...
        lightsOffMoment/60, ...
        lightsOnMoment/60, ...
        sleepOnsetTimepointafterscoringstart/60, ...
        sleepOffsetTimepointafterscoringstart/60, ...
        N1_bouts, ...
        N2_bouts, ...
        N3_bouts, ...
        S4_bouts, ...
        REM_bouts, ...
        Wake_bouts, ...
        MovementTime_bouts, ...
        Artifact_bouts, ...
        Unknown_bouts, ...
        SWS_bouts, ...
        NR_bouts, ...
        NRwithN1_bouts, ...
        S_bouts, ...
        nonS_bouts, ...
        WakeMT_bouts, ...
        NtransitionsTST, ...
        NtransitionsTST/(totalSleepPeriod/60), ...
        nArousals, ...
        nArousals_sleep_period, ...
        nArousalsByStage(1), ...
        nArousalsByStage(2), ...
        nArousalsByStage(3), ...
        nArousalsByStage(4), ...
        nArousalsByStage(5), ...
        nArousalsByStage(6), ...
        nArousals_sleep_period/(totalSleepPeriod/60), ...
        nArousals_sleep_period/(totalSleepInSleepPeriod/60), ...
        nArousalsByStage(1)/(N1Time/60), ...
        nArousalsByStage(2)/(N2Time/60), ...
        nArousalsByStage(3)/(N3Time/60), ...
        nArousalsByStage(4)/(S4Time/60), ...
        nArousalsByStage(5)/(REMtime/60), ...
        nArousalsByStage(6)/(WakeTime/60), ...
        overall_Duration_min_by_IntType_and_EpType(1,1,1), ...
        mean_Duration_min_by_IntType_and_EpType(1,1,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,1,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,1,1)/timewindow_min_by_IntType_and_EpType(1,1,2), ...
        overall_Duration_min_by_IntType_and_EpType(1,2,1), ...
        mean_Duration_min_by_IntType_and_EpType(1,2,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,2,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,2,1)/timewindow_min_by_IntType_and_EpType(1,2,2), ...
        overall_Duration_min_by_IntType_and_EpType(1,3,1), ...
        mean_Duration_min_by_IntType_and_EpType(1,3,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,3,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,3,1)/timewindow_min_by_IntType_and_EpType(1,3,2), ...
        overall_Duration_min_by_IntType_and_EpType(2,1,1), ...
        mean_Duration_min_by_IntType_and_EpType(2,1,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(2,1,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(2,1,1)/timewindow_min_by_IntType_and_EpType(2,1,2), ...
        overall_Duration_min_by_IntType_and_EpType(2,3,1), ...
        mean_Duration_min_by_IntType_and_EpType(2,3,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(2,3,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(2,3,1)/timewindow_min_by_IntType_and_EpType(2,3,2), ...
        overall_Duration_min_by_IntType_and_EpType(3,2,1), ...
        mean_Duration_min_by_IntType_and_EpType(3,2,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(3,2,1), ...
        Interruptions_in_bouts_by_IntType_and_EpType(3,2,1)/timewindow_min_by_IntType_and_EpType(3,2,2), ...
        overall_Duration_min_by_IntType_and_EpType(1,1,2), ...
        mean_Duration_min_by_IntType_and_EpType(1,1,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,1,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,1,2)/timewindow_min_by_IntType_and_EpType(1,1,2), ...
        overall_Duration_min_by_IntType_and_EpType(1,2,2), ...
        mean_Duration_min_by_IntType_and_EpType(1,2,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,2,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,2,2)/timewindow_min_by_IntType_and_EpType(1,2,2), ...
        overall_Duration_min_by_IntType_and_EpType(1,3,2), ...
        mean_Duration_min_by_IntType_and_EpType(1,3,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,3,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(1,3,2)/timewindow_min_by_IntType_and_EpType(1,3,2), ...
        overall_Duration_min_by_IntType_and_EpType(2,1,2), ...
        mean_Duration_min_by_IntType_and_EpType(2,1,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(2,1,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(2,1,2)/timewindow_min_by_IntType_and_EpType(2,1,2), ...
        overall_Duration_min_by_IntType_and_EpType(2,3,2), ...
        mean_Duration_min_by_IntType_and_EpType(2,3,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(2,3,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(2,3,2)/timewindow_min_by_IntType_and_EpType(2,3,2), ...
        overall_Duration_min_by_IntType_and_EpType(3,2,2), ...
        mean_Duration_min_by_IntType_and_EpType(3,2,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(3,2,2), ...
        Interruptions_in_bouts_by_IntType_and_EpType(3,2,2)/timewindow_min_by_IntType_and_EpType(3,2,2), ...
        60*Interruptions_in_bouts_by_IntType_and_EpType(3,2,2)/timewindow_min_by_IntType_and_EpType(3,2,2), ...
        timewindow_min_by_IntType_and_EpType(1,1,1), ...
        timewindow_min_by_IntType_and_EpType(1,2,1), ...
        timewindow_min_by_IntType_and_EpType(1,3,1), ...
        cycle_labels(iScoringCycle), ...
        cycle_complete(iScoringCycle), ...
        'VariableNames',{'scoringfile','standard','sleep_onset_definition','epoch_length_seconds',...
        'data_offset_min',...
        'scoring_duration_min',...
        'sleep_opportunity_on_to_off_min',...
        'lightsoff_to_lightson_min',...
        'total_sleep_period_duration_min','total_sleep_duration_in_sleep_period_min','sleep_eff_total_sleep_dur_in_sleep_prd_perc_of_sleep_prd','sleep_eff_total_sleep_dur_in_sleep_prd_perc_of_sleep_opport','sleep_onset_delay_min','N1_delay_min','N2_delay_min','SWS_delay_min','S4_delay_min','R_delay_min',...
        'N1_of_sleep_period_min','N2_of_sleep_period_min','N3_of_sleep_period_min','S4_of_sleep_period_min','R_of_sleep_period_min','Wake_after_sleep_onset_of_sleep_period_min','Movement_Time_of_sleep_period_min','Artifact_of_sleep_period_min','Unknown_of_sleep_period_min','SWS_of_sleep_period_min','NR_with_N1_of_sleep_period_min','NR_without_N1_of_sleep_period_min',...
        'N1_all_scoring_min','N2_all_scoring_min','N3_all_scoring_min','S4_all_scoring_min','R_all_scoring_min','Wake_all_scoring_min','Movement_Time_all_scoring_min','Artifact_all_scoring_min','Unknown_all_scoring_min','SWS_all_scoring_min','NR_with_N1_all_scoring_min','NR_without_N1_all_scoring_min',...
        'N1_perc_of_sleep_period','N2_perc_of_sleep_period','N3_perc_of_sleep_period','S4_perc_of_sleep_period','R_perc_of_sleep_period','Wake_after_sleep_onset_perc_of_sleep_period','Movement_Time_perc_of_sleep_period','Artifact_perc_of_sleep_period','Unknown_perc_of_sleep_period','SWS_perc_of_sleep_period','NR_with_N1_perc_of_sleep_period','NR_without_N1_perc_of_sleep_period',...
        'N1_perc_of_slp_prd_sleep_duration','N2_perc_of_slp_prd_sleep_duration','N3_perc_of_slp_prd_sleep_duration','S4_perc_of_slp_prd_sleep_duration','R_perc_of_slp_prd_sleep_duration','Wake_after_sleep_onset_perc_of_slp_prd_sleep_duration','Movement_Time_perc_of_slp_prd_sleep_duration','Artifact_perc_of_slp_prd_sleep_duration','Unknown_perc_of_slp_prd_sleep_duration','SWS_perc_of_slp_prd_sleep_duration','NR_with_N1_perc_of_slp_prd_sleep_duration','NR_without_N1_perc_of_slp_prd_sleep_duration',...
        'N1_without_excluded_min','N2_without_excluded_min','N3_without_excluded_min','S4_without_excluded_min','R_without_excluded_min','Wake_after_sleep_onset_without_excluded_min','Movement_Time_without_excluded_min', 'Artifact_without_excluded_min', 'Unknown_without_excluded_min', 'SWS_without_excluded_min', 'NR_with_N1_without_excluded_min','NR_without_N1_without_excluded_min',...
        'N1_without_excluded_perc_of_sleep_period','N2_without_excluded_perc_of_sleep_period','N3_without_excluded_perc_of_sleep_period','S4_without_excluded_perc_of_sleep_period','R_without_excluded_perc_of_sleep_period','Wake_after_sleep_onset_without_excluded_perc_of_sleep_period','Movement_Time_without_excluded_perc_of_sleep_period', 'Artifact_without_excluded_perc_of_sleep_period', 'Unknown_without_excluded_perc_of_sleep_period','SWS_without_excluded_perc_of_sleep_period','NR_with_N1_without_excluded_perc_of_sleep_period','NR_without_N1_without_excluded_perc_of_sleep_period',...
        'N1_without_excluded_perc_of_slp_prd_sleep_duration','N2_without_excluded_perc_of_slp_prd_sleep_duration','N3_without_excluded_perc_of_slp_prd_sleep_duration','S4_without_excluded_perc_of_slp_prd_sleep_duration','R_without_excluded_perc_of_slp_prd_sleep_duration','Wake_after_SO_without_excluded_perc_of_slp_prd_slp_duration','Movement_Time_without_excluded_perc_of_slp_prd_sleep_duration', 'Artifact_without_excluded_perc_of_slp_prd_sleep_duration', 'Unknown_without_excluded_perc_of_slp_prd_sleep_duration','SWS_without_excluded_perc_of_slp_prd_sleep_duration','NR_with_N1_without_excluded_perc_of_slp_prd_sleep_duration','NR_without_N1_without_excluded_perc_of_slp_prd_sleep_duration',...
        'N1_before_sleep_onset_min','N2_before_sleep_onset_min','N3_before_sleep_onset_min','S4_before_sleep_onset_min','R_before_sleep_onset_min','Wake_before_sleep_onset_min','Movement_Time_before_sleep_onset_min', 'Artifact_before_sleep_onset_min', 'Unknown_before_sleep_onset_min','SWS_before_sleep_onset_min','NR_before_sleep_onset_without_N1_min',...
        'longest_WASO_period_of_sleep_period_min','longest_WASO_period_after_latency_after_SlpOn_min',...
        'allowed_sleep_onset_before_sleep_opportunity',...
        'allowed_sleep_offset_after_sleep_opportunity',...
        'sleep_opportunity_on_exists',...
        'sleep_opportunity_off_exists',...
        'lights_off_exists',...
        'lights_on_exists',...
        'sleep_opportunity_on_min',...
        'sleep_opportunity_off_min',...
        'lights_off_min',...
        'lights_on_min',...
        'sleep_onset_time_after_scoring_start_min',...
        'sleep_offset_time_after_scoring_start_min',...
        'N1_number_of_bouts_in_sleep_period', ...
        'N2_number_of_bouts_in_sleep_period', ...
        'N3_number_of_bouts_in_sleep_period', ...
        'S4_number_of_bouts_in_sleep_period', ...
        'REM_number_of_bouts_in_sleep_period', ...
        'Wake_number_of_bouts_in_sleep_period', ...
        'MovementTime_number_of_bouts_in_sleep_period', ...
        'Artifact_number_of_bouts_in_sleep_period', ...
        'Unknown_number_of_bouts_in_sleep_period', ...
        'SWS_number_of_bouts_in_sleep_period', ...
        'NRwithN1_number_of_bouts_in_sleep_period', ...
        'NR_number_of_bouts_in_sleep_period', ...
        'sleep_period_number_of_bouts', ...
        'non_sleep_period_number_of_bouts', ...
        'WakeAndMovementTime_number_of_bouts_in_sleep_period', ...
        'stage_transitions_in_sleep_period', ...
        'stage_transitions_density_per_min_in_sleep_period', ...
        'arousals_number', ...
        'arousals_number_in_sleep_period', ...
        'arousals_number_N1_in_sleep_period', ...
        'arousals_number_N2_in_sleep_period', ...
        'arousals_number_N3_in_sleep_period', ...
        'arousals_number_S4_in_sleep_period', ...
        'arousals_number_R_in_sleep_period',...
        'arousals_number_WASO_in_sleep_period', ...
        'arousals_density_per_min_sleep_period', ...
        'arousals_density_per_min_sleep_in_sleep_period', ...
        'arousals_density_per_min_N1_in_sleep_period', ...
        'arousals_density_per_min_N2_in_sleep_period', ...
        'arousals_density_per_min_N3_in_sleep_period', ...
        'arousals_density_per_min_S4_in_sleep_period', ...
        'arousals_density_per_min_R_in_sleep_period',...
        'arousals_density_per_min_WASO_in_sleep_period', ...
        'interrupts_W_MT_overall_durations_min_in_NR_epsds', ...
        'interrupts_W_MT_mean_durations_min_in_NR_epsds', ...
        'interrupts_W_MT_number_in_NR_epsds', ...
        'interrupts_W_MT_density_per_min_in_NR_epsds', ...
        'interrupts_W_MT_overall_durations_min_in_R_epsds', ...
        'interrupts_W_MT_mean_durations_min_in_R_epsds', ...
        'interrupts_W_MT_number_in_R_epsds', ...
        'interrupts_W_MT_density_per_min_in_R_epsds', ...
        'interrupts_W_MT_overall_durations_min_in_ascNR_epsds', ...
        'interrupts_W_MT_mean_durations_min_in_ascNR_epsds', ...
        'interrupts_W_MT_number_in_ascNR_epsds', ...
        'interrupts_W_MT_density_per_min_in_ascNR_epsds', ...
        'interrupts_W_MT_N1_overall_durations_min_in_NR_epsds', ...
        'interrupts_W_MT_N1_mean_durations_min_in_NR_epsds', ...
        'interrupts_W_MT_N1_number_in_NR_epsds', ...
        'interrupts_W_MT_N1_density_per_min_in_NR_epsds', ....
        'interrupts_W_MT_N1_overall_durations_min_in_ascNR_epsds', ...
        'interrupts_W_MT_N1_mean_durations_min_in_ascNR_epsds', ...
        'interrupts_W_MT_N1_number_in_ascNR_epsds', ...
        'interrupts_W_MT_N1_density_per_min_in_ascNR_epsds', ...
        'interrupts_W_MT_NR_overall_durations_min_in_R_epsds', ...
        'interrupts_W_MT_NR_mean_durations_min_in_R_epsds', ...
        'interrupts_W_MT_NR_number_in_R_epsds', ...
        'interrupts_W_MT_NR_density_per_min_in_R_epsds', ...
        'interrupts_W_MT_arsls_overall_durations_min_in_NR_epsds', ...
        'interrupts_W_MT_arsls_mean_durations_min_in_NR_epsds', ...
        'interrupts_W_MT_arsls_number_in_NR_epsds', ...
        'interrupts_W_MT_arsls_density_per_min_in_NR_epsds', ...
        'interrupts_W_MT_arsls_overall_durations_min_in_R_epsds', ...
        'interrupts_W_MT_arsls_mean_durations_min_in_R_epsds', ...
        'interrupts_W_MT_arsls_number_in_R_epsds', ...
        'interrupts_W_MT_arsls_density_per_min_in_R_epsds', ...
        'interrupts_W_MT_arsls_overall_durations_min_in_ascNR_epsds', ...
        'interrupts_W_MT_arsls_mean_durations_min_in_ascNR_epsds', ...
        'interrupts_W_MT_arsls_number_in_ascNR_epsds', ...
        'interrupts_W_MT_arsls_density_per_min_in_ascNR_epsds', ...
        'interrupts_W_MT_N1_arsls_overall_durations_min_in_NR_epsds', ...
        'interrupts_W_MT_N1_arsls_mean_durations_min_in_NR_epsds', ...
        'interrupts_W_MT_N1_arsls_number_in_NR_epsds', ...
        'interrupts_W_MT_N1_arsls_density_per_min_in_NR_epsds', ....
        'interrupts_W_MT_N1_arsls_overall_durations_min_in_ascNR_epsds', ...
        'interrupts_W_MT_N1_arsls_mean_durations_min_in_ascNR_epsds', ...
        'interrupts_W_MT_N1_arsls_number_in_ascNR_epsds', ...
        'interrupts_W_MT_N1_arsls_density_per_min_in_ascNR_epsds', ...
        'interrupts_W_MT_NR_arsls_overall_durations_min_in_R_epsds', ...
        'interrupts_W_MT_NR_arsls_mean_durations_min_in_R_epsds', ...
        'interrupts_W_MT_NR_arsls_number_in_R_epsds', ...
        'interrupts_W_MT_NR_arsls_density_per_min_in_R_epsds', ...
        'REM_fragmentation_density_per_hour_in_R_epsds', ...
        'cycle_overall_duration_min_NR_epsds', ...
        'cycle_overall_duration_min_R_epsds', ...
        'cycle_overall_duration_min_ascNR_epsds', ...
        'cycle_label',...
        'cycle_complete'}...
        );

    if iScoringCycle == 1
        res.table = table_tmp;
    else
        res.table = cat(1,res.table,table_tmp);
    end

end
fprintf([functionname ' function finished\n']);

toc
memtoc
end
