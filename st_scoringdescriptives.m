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
%   cfg.cycle          = string or integer value from 1 to Inf. 
%                        If set, then the descriptives is derived from 
%                        only this cycle.
%                        possible values are:
%                          'presleep' = prior to first cycle
%                          1 = first cycle
%                          2 = second cycle
%                          ... 
%                          'last' = last complete cycle
%                          'postsleep' = post last complete cycle
%                          'all' = post last complete cycle
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
    
    indcycle = cfg.cycle;
    switch cfg.cycle
        case 'last'
            indcycle = strcmp(cycle_labels,'last');
        case 'presleep'
            indcycle = strcmp(cycle_labels,'presleep');
        case 'postsleep'
            indcycle = strcmp(cycle_labels,'postsleep');
        case 'complete'
            indcycle = cycle_complete;
        case 'all'
            indcycle = 1:numel(scoring_cycles);
    end
    scoring_cycles = scoring_cycles(indcycle);
    cycle_complete = cycle_complete(indcycle);
    cycle_labels = cycle_labels(indcycle);
    scoring_duration_min = scoring_duration_min(indcycle);
    
    data_offset_min_cycle = nan(1,numel(cfg2.end));
    for iSc = 1:numel(cfg2.end)
        if ~isempty(scoring_duration_min)
            data_offset_min_cycle(iSc) = data_offset_min + sum(scoring_duration_min(1:(iSc-1)));
        end
    end
    data_offset_min = data_offset_min_cycle(indcycle);
end


res = [];
res.ori = getfunctionname();
res.type = 'descriptive';
res.cfg = cfg;

for iScoringCycle = 1:numel(scoring_cycles)
    scoring = scoring_cycles{iScoringCycle};
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
    stagesCombo = {{'N1','S1'}, {'N2','S2'}, {'N3','S3'}, {'S4'}, {'R'}, {'W'}};
    nArousalsByStage = nan(size(stagesCombo));

    fprintf([functionname ' function initialized\n']);
    
    
    dummySampleRate = 100;
    %lightsOffSample = round(lightsOffMoment*dummySampleRate);
    epochLengthSamples = scoring.epochlength * dummySampleRate;
    
    epochs = cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0);
    hypnStages = [cellfun(@sleepStage2str,epochs,'UniformOutput',0) ...
        cellfun(@sleepStage2str_alt,epochs,'UniformOutput',0) ...
        cellfun(@sleepStage2str_alt2,epochs,'UniformOutput',0) ...
        cellfun(@sleepStage2str_alt3,epochs,'UniformOutput',0)];
    
    
    hypnEpochs = 1:numel(epochs);
    hypnEpochsBeginsSamples = (((hypnEpochs - 1) * epochLengthSamples) + 1)';
    hypnEpochsEndsSamples   = (hypnEpochs * epochLengthSamples)';
    
    [onsetnumber lastsleepstagenumber onsetepoch lastsleepstage allowedsleeponsetbeforesleepopon] = st_sleeponset(cfg,scoring);
    
    if isempty(onsetepoch)
        onsetnumber = NaN;
    end
    
    hasLightsOn = false;
    lightsOnMoment = lastsleepstagenumber*scoring.epochlength;
    lightsOnToLightsOff = NaN;
    if isfield(scoring, 'lightson')
        if ~isnan(scoring.lightson)
            hasLightsOn = true;
            lightsOnMoment = scoring.lightson;
            lightsOnToLightsOff = (lightsOnMoment-lightsOffMoment);
        else
            ft_warning('The lights on moment was NaN in the scoring structure.\n The end of sleep is thus assumed as lights on.');
        end
    else
        ft_warning('The lights on moment was not provided in the scoring structure.\n The end of sleep is thus assumed as lights on.');
    end
    
    
    hasSleepOpportunityOff = false;
    sleepOpportunityOffMoment = lastsleepstagenumber*scoring.epochlength;
    sleepOpportunityOnTosleepOpportunityOff = NaN;
    if isfield(scoring, 'sleepopoff')
        if ~isnan(scoring.sleepopoff)
            hasSleepOpportunityOff = true;
            sleepOpportunityOffMoment = scoring.sleepopoff;
            sleepOpportunityOnTosleepOpportunityOff = (sleepOpportunityOffMoment-sleepOpportunityOnMoment);
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
    
    if allowedsleeponsetbeforesleepopon
        ft_warning('sleep onset was allowed to be before sleep opportunity moment and thus sleep onset latency will be NaN.')
        sleepOnsetTime = NaN;
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
        hypnTST_excluded = scoring.excluded(onsetnumber:lastsleepstagenumber);
        
        hypnStagesTST = hypnStages(onsetnumber:lastsleepstagenumber,:);
        hypnStagesPreSleepOnset = hypnStages(1:preOnsetCandidate,:);
        hypnEpochsTST = hypnEpochs(onsetnumber:lastsleepstagenumber);
        hypnEpochsBeginsSamplesTST = hypnEpochsBeginsSamples(onsetnumber:lastsleepstagenumber,:);
        hypnEpochsEndsSamplesTST = hypnEpochsEndsSamples(onsetnumber:lastsleepstagenumber,:);
        
        totalSleepPeriod = (length(onsetnumber:lastsleepstagenumber))*scoring.epochlength;
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
    
    
    if isfield(scoring,'arousals')
        nArousals = size(scoring.arousals,1);
        
        for iStages = 1:numel(stagesCombo)
            

            if ~isnan(onsetnumber)
                cfg_cut = [];
                cfg_cut.start = onsetnumber;
                cfg_cut.end = lastsleepstagenumber;
                scorings_sleepperiod = st_cutscoring(cfg_cut,scoring);
                scorings_sleepperiod = scorings_sleepperiod{1};
            else
                scorings_sleepperiod = scoring;
                scorings_sleepperiod.epochs = scorings_sleepperiod.epochs(1:0);
                scorings_sleepperiod.excluded = scorings_sleepperiod.excluded(1:0);
            end
            
            cfg = [];
            cfg.stages = stagesCombo{iStages};
            cfg.considerexcluded = 'no';
            [begins_epoch ends_epoch] = st_select_scoring(cfg,scorings_sleepperiod);
            if ~isempty(begins_epoch)
                cfg.start = begins_epoch;
                cfg.end = ends_epoch;
                scorings_cut = st_cutscoring(cfg,scorings_sleepperiod);
                
                scorings_cut_appended = st_append_scoring(scorings_cut{:});
                nArousalsByStage(iStages) = size(scorings_cut_appended.arousals,1);
            end
        end
    end
    
    hypnEpochs = 1:numel(scoring.epochs);
    
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'N1'})),1);
    N1_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'N2'})),1);
    N2_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'N3'})),1);
    N3_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'S4'})),1);
    S4_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'R'})),1);
    REM_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'W'})),1);
    Wake_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'MT'})),1);
    MovementTime_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'A'})),1);
    Artifact_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'?'})),1);
    Unknown_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'N3','S4'})),1);
    SWS_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'N2','N3','S4'})),1);
    NR_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,1),{'N1','N2','N3','S4'})),1);
    NRwithN1_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,4),{'S'})),1);
    S_bouts = numel(consecBegins);
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(hypnEpochs(ismember(hypnStagesTST(:,4),{'W'})),1);
    nonS_bouts = numel(consecBegins);
    
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
    
    table_tmp = table(...
        {scoringfile}, ...
        {scoring.standard}, ...
        scoring.epochlength, ...
        totalSleepPeriod/60, ...
        (totalSleepPeriod-WakeTime)/60, ...
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
        100*N1Time/(totalSleepPeriod-WakeTime), ...
        100*N2Time/(totalSleepPeriod-WakeTime), ...
        100*N3Time/(totalSleepPeriod-WakeTime), ...
        100*S4Time/(totalSleepPeriod-WakeTime), ...
        100*REMtime/(totalSleepPeriod-WakeTime), ...
        100*WakeTime/(totalSleepPeriod-WakeTime), ...
        100*MovementTime/(totalSleepPeriod-WakeTime), ...
        100*ArtifactTime/(totalSleepPeriod-WakeTime), ...
        100*UnknownTime/(totalSleepPeriod-WakeTime), ...
        100*SWStime/(totalSleepPeriod-WakeTime), ...
        100*(NRtime+N1Time)/(totalSleepPeriod-WakeTime), ...
        100*NRtime/(totalSleepPeriod-WakeTime), ...
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
        100*N1Time_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*N2Time_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*N3Time_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*S4Time_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*REMtime_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*WakeTime_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*MovementTime_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*ArtifactTime_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*UnknownTime_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*SWStime_without_excluded/(totalSleepPeriod-WakeTime), ...
        100*(NRtime_without_excluded+N1Time_without_excluded)/(totalSleepPeriod-WakeTime), ...
        100*NRtime_without_excluded/(totalSleepPeriod-WakeTime), ...
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
        lightsOnToLightsOff/60, ...
        sleepOpportunityOnTosleepOpportunityOff/60, ...
        allowedsleeponsetbeforesleepopon, ...
        hasSleepOpportunityOn, ...
        hasSleepOpportunityOff, ...
        hasLightsOff, ...
        hasLightsOn, ...
        sleepOpportunityOnMoment/60, ...
        sleepOpportunityOffMoment/60, ...
        lightsOffMoment/60, ...
        lightsOnMoment/60, ...
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
        nArousals, ...
        nArousalsByStage(1), ...
        nArousalsByStage(2), ...
        nArousalsByStage(3), ...
        nArousalsByStage(4), ...
        nArousalsByStage(5), ...
        nArousalsByStage(6), ...
        nArousalsByStage(1)/(N1Time/60), ...
        nArousalsByStage(2)/(N2Time/60), ...
        nArousalsByStage(3)/(N3Time/60), ...
        nArousalsByStage(4)/(S4Time/60), ...
        nArousalsByStage(5)/(REMtime/60), ...
        nArousalsByStage(6)/(WakeTime/60), ...
        cycle_labels(iScoringCycle), ...
        cycle_complete(iScoringCycle), ...
        scoring_duration_min(iScoringCycle),...
        data_offset_min(iScoringCycle),...
        'VariableNames',{'scoringfile','standard','epoch_length_seconds','total_sleep_period_duration_min','total_sleep_duration_in_sleep_period_min','sleep_onset_min','N1_delay_min','N2_delay_min','SWS_delay_min','S4_delay_min','R_delay_min'...
        ,'N1_min','N2_min','N3_min','S4_min','R_min','Wake_after_sleep_onset_min','Movement_Time_min','Artifact_min','Unknown_min','SWS_min','NR_with_N1_min','NR_without_N1_min'...
        ,'N1_perc_of_sleep_period','N2_perc_of_sleep_period','N3_perc_of_sleep_period','S4_perc_of_sleep_period','R_perc_of_sleep_period','Wake_after_sleep_onset_perc_of_sleep_period','Movement_Time_perc_of_sleep_period','Artifact_perc_of_sleep_period','Unknown_perc_of_sleep_period','SWS_perc_of_sleep_period','NR_with_N1_perc_of_sleep_period','NR_without_N1_perc_of_sleep_period'...
        ,'N1_perc_of_sleep_duration','N2_perc_of_sleep_duration','N3_perc_of_sleep_duration','S4_perc_of_sleep_duration','R_perc_of_sleep_duration','Wake_after_sleep_onset_perc_of_sleep_duration','Movement_Time_perc_of_sleep_duration','Artifact_perc_of_sleep_duration','Unknown_perc_of_sleep_duration','SWS_perc_of_sleep_duration','NR_with_N1_perc_of_sleep_duration','NR_without_N1_perc_of_sleep_duration'...
        ,'N1_without_excluded_min','N2_without_excluded_min','N3_without_excluded_min','S4_without_excluded_min','R_without_excluded_min','Wake_after_sleep_onset_without_excluded_min','Movement_Time_without_excluded_min', 'Artifact_without_excluded_min', 'Unknown_without_excluded_min', 'SWS_without_excluded_min', 'NR_with_N1_without_excluded_min','NR_without_N1_without_excluded_min'...
        ,'N1_without_excluded_perc_of_sleep_period','N2_without_excluded_perc_of_sleep_period','N3_without_excluded_perc_of_sleep_period','S4_without_excluded_perc_of_sleep_period','R_without_excluded_perc_of_sleep_period','Wake_after_sleep_onset_without_excluded_perc_of_sleep_period','Movement_Time_without_excluded_perc_of_sleep_period', 'Artifact_without_excluded_perc_of_sleep_period', 'Unknown_without_excluded_perc_of_sleep_period','SWS_without_excluded_perc_of_sleep_period','NR_with_N1_without_excluded_perc_of_sleep_period','NR_without_N1_without_excluded_perc_of_sleep_period'...
        ,'N1_without_excluded_perc_of_sleep_duration','N2_without_excluded_perc_of_sleep_duration','N3_without_excluded_perc_of_sleep_duration','S4_without_excluded_perc_of_sleep_duration','R_without_excluded_perc_of_sleep_duration','Wake_after_sleep_onset_without_excluded_perc_of_sleep_duration','Movement_Time_without_excluded_perc_of_sleep_duration', 'Artifact_without_excluded_perc_of_sleep_duration', 'Unknown_without_excluded_perc_of_sleep_duration','SWS_without_excluded_perc_of_sleep_duration','NR_with_N1_without_excluded_perc_of_sleep_duration','NR_without_N1_without_excluded_perc_of_sleep_duration'...
        ,'N1_before_sleep_onset_min','N2_before_sleep_onset_min','N3_before_sleep_onset_min','S4_before_sleep_onset_min','R_before_sleep_onset_min','Wake_before_sleep_onset_min','Movement_Time_before_sleep_onset_min', 'Artifact_before_sleep_onset_min', 'Unknown_before_sleep_onset_min','SWS_before_sleep_onset_min','NR_before_sleep_onset_without_N1_min'...
        ,'longest_WASO_period_min','longest_WASO_period_after_latency_after_so_min',...
        'lightsoff_to_lightson_min',...
        'sleep_opportunity_on_to_off_min',...
        'allowed_sleep_onset_before_sleep_opportunity',...
        'sleep_opportunity_on_present',...
        'sleep_opportunity_off_present',...
        'lights_off_present',...
        'lights_on_present',...
        'sleep_opportunity_on',...
        'sleep_opportunity_off',...
        'lights_off',...
        'lights_on',...
        'N1_number_of_bouts', ...
        'N2_number_of_bouts', ...
        'N3_number_of_bouts', ...
        'S4_number_of_bouts', ...
        'REM_number_of_bouts', ...
        'Wake_number_of_bouts', ...
        'MovementTime_number_of_bouts', ...
        'Artifact_number_of_bouts', ...
        'Unknown_number_of_bouts', ...
        'SWS_number_of_bouts', ...
        'NRwithN1_number_of_bouts', ...
        'NR_number_of_bouts',...
        'sleep_period_number_of_bouts',...
        'non_sleep_period_number_of_bouts',...
        'arousals_number', ...
        'arousals_number_N1', ...
        'arousals_number_N2', ...
        'arousals_number_N3', ...
        'arousals_number_S4', ...
        'arousals_number_R',...
        'arousals_number_WASO', ...
        'arousals_density_per_min_N1', ...
        'arousals_density_per_min_N2', ...
        'arousals_density_per_min_N3', ...
        'arousals_density_per_min_S4', ...
        'arousals_density_per_min_R',...
        'arousals_density_per_min_WASO', ...
        'cycle_label',...
        'cycle_complete',...
        'scoring_duration_min',...
        'data_offset_min'}...
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
