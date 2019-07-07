function [startepochs endepochs Rstarts Rends otherstarts otherends] = st_sleepcycles(cfg,scoring)

% ST_SLEEPCYCLES find the sleep cycles
%
% Use as
%   [startepochs endepochs] = st_sleepcycles(cfg,scoring)
%   [startepochs endepochs Rstarts Rends otherstarts otherends] = st_sleepcycles(cfg,scoring)
%
% Configutation parameter can be empty, e.g. cfg = []
%
% Optional configuration parameters are
%   cfg.smoothepochs   = number of epochs to smooth and find REM, default
%                        is the amount of epochs that fill 15 minutes
%   cfg.sleeponsetdef  = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                        'NR' or 'NR' or 'XR', see below for details (default = 'N1_XR')
%
% See also ST_READ_SCORING



cfg.smoothepochs       = ft_getopt(cfg, 'smoothepochs', 15*60/scoring.epochlength);
cfg.sleeponsetdef      = ft_getopt(cfg, 'sleeponsetdef', 'N1_XR');

hasLightsOff = false;
if isfield(scoring, 'lightsoff')
    hasLightsOff = true;
end

nEpochs = numel(scoring.epochs);

hypnStages = [cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0) ...
    cellfun(@sleepStage2str_alt,scoring.epochs','UniformOutput',0) ...
    cellfun(@sleepStage2str_alt2,scoring.epochs','UniformOutput',0)];

[onsetCandidateIndex, preOffsetCandidate, onsetepoch] = st_sleeponset(cfg,scoring);
if isempty(onsetepoch)
    onsetCandidateIndex = nEpochs;
end

if isempty(preOffsetCandidate)
    preOffsetCandidate = nEpochs;
end

epochs = hypnStages(:,3)';

[Rstarts, Rends] = getCycleStartEndsByLabel('R',cfg.smoothepochs,epochs);

if isempty(Rends)
    Rends = preOffsetCandidate;
end

cycleStarts = onsetCandidateIndex;
cycleEnds = Rends(1);

for iREM_ends = 2:numel(Rends)
    cycleStarts = [cycleStarts; Rends(iREM_ends-1)+1];
    cycleEnds   = [cycleEnds;   Rends(iREM_ends)];
end

startepochs = cycleStarts;
endepochs = cycleEnds;


otherstarts = cycleStarts;
otherends = Rstarts-1;

if Rends(end) < preOffsetCandidate
    otherstarts = [otherstarts; Rends(end)+1];
    otherends   = [otherends; preOffsetCandidate];
end

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



