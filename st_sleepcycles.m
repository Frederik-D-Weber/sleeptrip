function [res_cycle] = st_sleepcycles(cfg,scoring)

% ST_SLEEPCYCLES find the sleep cycles
%
% Use as
%   res_cycle = st_sleepcycles(cfg, scoring)
%
% Configutation parameter can be empty, e.g. cfg = []
%
% Optional configuration parameters are
%   either
%   cfg.smoothepochs   = number of epochs to smooth and find REM, default
%                        is the amount of epochs that fill 15 minutes
%   or
%   cfg.smoothminutes  = minutes to smooth and find REM, default
%                        is 15 minutes and. This option overwrites
%                        cfg.smoothepochs if present
%
%   cfg.sleeponsetdef  = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                        'NR' or 'N2R' or 'XR' or 'AASM' or 'X2R' or 
%                        'N2' or 'N3' or 'SWS' or 'S4' or 'R, see ST_SLEEPONSET for details (default = 'AASM')
%   cfg.maxduration    = duration in minutes a sleep cycle can maximally
%                        have cannot be lower than 2. (default = Inf) 
%   cfg.completesleepbuffer  = duration of sleep in minutes a 
%                              sleep cycle has to be followed to be
%                              considered complete (default = 5)
%
% result table has columns with: cycle startepochs endepochs Rstarts Rends NRstarts NRends
%
% See also ST_READ_SCORING, ST_SLEEPONSET

timerVal = tic;
memtic
st = dbstack;
functionname = st.name;

fprintf([functionname ' function started\n']);


cfg.smoothepochs      = ft_getopt(cfg, 'smoothepochs', 15*60/scoring.epochlength);

if isfield(cfg, 'smoothminutes')
    cfg.smoothepochs       = ft_getopt(cfg, 'smoothepochs', floor(cfg.smoothminutes*60/scoring.epochlength));
else
    cfg.smoothepochs       = ft_getopt(cfg, 'smoothepochs', floor(15*60/scoring.epochlength));
end

cfg.maxduration        = ft_getopt(cfg, 'maxduration', Inf);
cfg.completesleepbuffer        = ft_getopt(cfg, 'completesleepbuffer', 5);

cfg.sleeponsetdef      = upper(ft_getopt(cfg, 'sleeponsetdef', 'AASM'));

maxDurationEpochs = max(2,floor(cfg.maxduration*60/scoring.epochlength));
completesleepbufferDurationEpochs = max(2,floor(cfg.completesleepbuffer*60/scoring.epochlength));


hasLightsOff = false;
if isfield(scoring, 'lightsoff')
    hasLightsOff = true;
end

nEpochs = numel(scoring.epochs);

hasSleep = true;
[onsetCandidateIndex, preOffsetCandidate, onsetepoch] = st_sleeponset(cfg,scoring);
if isempty(onsetepoch) % no sleep onset found
    hasSleep = false;
    onsetCandidateIndex = nEpochs;
end

if isempty(preOffsetCandidate)
    preOffsetCandidate = nEpochs;
end

% hypnStages = [cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0) ...
%     cellfun(@sleepStage2str_alt,scoring.epochs','UniformOutput',0) ...
%     cellfun(@sleepStage2str_alt2,scoring.epochs','UniformOutput',0) ...
%     cellfun(@sleepStage2str_alt3,scoring.epochs','UniformOutput',0)];
%
% epochs = hypnStages(:,3)';

epochs = cellfun(@sleepStage2str_alt2,scoring.epochs','UniformOutput',0)';

[Rstarts, Rends] = getCycleStartEndsByLabel('R',cfg.smoothepochs,epochs);

hasREM = true;
if isempty(Rstarts)
    hasREM = false;
    ft_warning('There is no REM at all!')
    
else
    
    idx = Rends > onsetCandidateIndex;
    Rends = Rends(idx);
    Rstarts = Rstarts(idx);
    
    if isempty(Rstarts)
        ft_warning('There is no REM after sleep onset!')
        hasREM = false;
    end

end
    
if isempty(Rends)
    hasREM = false;
    Rends = preOffsetCandidate;
end


if ~isempty(Rstarts)
    Rstarts(Rstarts < onsetCandidateIndex) = onsetCandidateIndex;

    if Rstarts(1) == onsetCandidateIndex
        ft_warning('First sleep cycle starts with REM!')
    end
end

cycleStarts = onsetCandidateIndex;
cycleEnds = Rends(1);

for iREM_ends = 2:numel(Rends)
    cycleStarts = [cycleStarts; Rends(iREM_ends-1)+1];
    cycleEnds   = [cycleEnds;   Rends(iREM_ends)];
end



cycleLengths = (cycleEnds - cycleStarts + 1);

whichCycleLengthsExcede = find(cycleLengths > maxDurationEpochs);
for iCycleIndex = 1:numel(whichCycleLengthsExcede)
    iCycle = whichCycleLengthsExcede(iCycleIndex);
    ft_warning(['Cyle ' num2str(iCycle) ' exceded the maximal duration of the epochs'])
    
%     %convert the sleep stages to hypnogram numbers
%     hypn = [cellfun(@(st) sleepStage2hypnNum(st,~istrue(cfg.plotunknown),false),scoring.epochs','UniformOutput',1), scoring.excluded'];
%     
%     
%     smoothingepochlength = cfg.smoothepochs;
%     if exist('smooth','file') == 2
%         smoothhyp = smooth(hypn(:,1),smoothingepochlength,'moving');
%     else
%         smoothhyp = smoothwd(hypn(:,1),smoothingepochlength)';
%     end
%     
%     hypn_part = smoothhyp(startepochs(iCycle):cycleEnds(iCycle));
%     
%     
    
    

end





startepochs = cycleStarts;
endepochs = cycleEnds;


NRstarts = cycleStarts;
if hasREM
    NRends = Rstarts-1;
else
    NRends = endepochs;
end

if Rends(end) < preOffsetCandidate
    NRstarts = [NRstarts; Rends(end)+1];
    NRends   = [NRends; preOffsetCandidate];
end


if ~hasREM
    Rstarts(1:end) = NaN;
    Rends(1:end) = NaN;
end

maxCycle = max([numel(startepochs),numel(endepochs),numel(Rstarts),numel(Rends),numel(NRstarts),numel(NRends)]);

startepochs((end+1):maxCycle) = NaN;
endepochs((end+1):maxCycle) = NaN;
Rstarts((end+1):maxCycle) = NaN;
Rends((end+1):maxCycle) = NaN;
NRstarts((end+1):maxCycle) = NaN;
NRends((end+1):maxCycle) = NaN;


startepochs = startepochs(:).'';
endepochs = endepochs(:).'';
Rstarts = Rstarts(:).'';
Rends = Rends(:).'';
NRstarts = NRstarts(:).'';
NRends = NRends(:).'';

res_cycle = [];
res_cycle.ori = functionname;
res_cycle.type = 'cycle';
res_cycle.cfg = cfg;
res_cycle.table = table(...
    (1:maxCycle)', startepochs(:), endepochs, NRstarts, NRends, Rstarts, Rends,...
    'VariableNames',{...
    'cycle', 'startepoch', 'endepoch', 'NRstartepoch', 'NRendepoch', 'Rstartepoch', 'Rendepoch'}...
    );

res_cycle.table.durationepochs = res_cycle.table.endepoch - res_cycle.table.startepoch + 1;
res_cycle.table.durationNRepochs = res_cycle.table.NRendepoch - res_cycle.table.NRstartepoch + 1;
res_cycle.table.durationRepochs = res_cycle.table.Rendepoch - res_cycle.table.Rstartepoch + 1;

res_cycle.table.NRstartepoch(res_cycle.table.durationNRepochs == 0) = NaN;
res_cycle.table.NRendepoch(res_cycle.table.durationNRepochs == 0) = NaN;

if ~hasSleep
    ft_warning('There was no sleep for sleep cycles to be detected, result has an empty table now.')
    res_cycle.table = res_cycle.table(1:0,:);
end

nCycles = size(res_cycle.table,1);
res_cycle.table.complete = logical(zeros(nCycles,1));
for iCycle = 1:(nCycles-1)
    if res_cycle.table.durationepochs(iCycle+1) >= completesleepbufferDurationEpochs;
        res_cycle.table.complete(iCycle) = logical(1);
    end
end

for iCycle = 1:(1)
iCycle
end

fprintf([functionname ' function finished\n']);
toc(timerVal)
memtoc
end

function [starts, ends] = getCycleStartEndsByLabel(label,max_merge_inbetween_stages,stages)
stages = getSmoothedLabels(label,max_merge_inbetween_stages,stages);
label_logical = strcmp(stages,label);
if label_logical(end) == 1
    label_logical = [label_logical 0];
end
onoff = diff(label_logical);

starts = find(onoff == 1)+1;
ends   = find(onoff == -1);
starts = starts';
ends = ends';
end

function stages = getSmoothedLabels(label,max_merge_inbetween_stages,stages)
if (max_merge_inbetween_stages > 0)
    is_stage = strcmp(stages,label);
    for iS = 2:length(is_stage)
        max_merge = min(iS-2,max_merge_inbetween_stages);
        if (max_merge_inbetween_stages < 1)
            continue
        end
        ind_merge_local = ((iS-max_merge-1):(iS-1));
        if (is_stage(iS) && any(is_stage(ind_merge_local)))
            first_prev_label_in_range = iS - length(ind_merge_local) - 1 + (find(is_stage(ind_merge_local),1));
            ind_change = first_prev_label_in_range:(iS-1);
            for iChange = ind_change
                stages{iChange} = label;
            end
        end
    end
end

end


% function [starts, ends] = startStop(stages,label)
% starts = [];
% ends = [];
% is_stage = strcmp(stages,label);
% if (length(is_stage) > 1)
%     run_len = seqle(is_stage);
%     
%     run_len.lengths.cumsum = cumsum(run_len$lengths)
%     run_select = which(run_len$values == T & run_len$lengths >= 1)
%     
%     ends = run_len.lengths.cumsum[run_select]
%     
%     newindex = ifelse(run_select>1, run_select-1, 0);
%     starts = run_len.lengths.cumsum[newindex] + 1;
%     if (0 %in% newindex) starts = c(1,starts)
%         
%     elseif (length(is_stage) == 1)
%         
%         if (is_stage(1))
%             starts = 1;
%             ends = 1;
%         end
%         
%     end
%     
% end
% end
% 
% 



