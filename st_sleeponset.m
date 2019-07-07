function [onsetnumber preoffsetnumber onsetepoch] = st_sleeponset(cfg,scoring)
% 
% ST_SLEEPONSET determines the sleep onset of a sleep scoring
% Use as
%   [onsetnumber, preoffsetnumber, onsetepoch] = st_sleeponset(cfg,scoring)
%   [onsetnumber, preoffsetnumber] = st_sleeponset(cfg,scoring)
%   [onsetnumber] = st_sleeponset(cfg,scoring)
%
% Configutation parameter can be empty, e.g. cfg = []
%
% Optional configuration parameters are
%   cfg.sleeponsetdef  = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                        'NR' or 'NR' or 'XR', see below for details (default = 'N1_XR')
%
%  Here are the possible sleep onset definitions, where NR referes to
%  non-REM, R to REM and XR to any non-REM or REM sleep stage.
%  N1     starting with first N1 is sleep onset, nothing else
%  N1_NR  starting with first N1 followed directly by either N2, N3 (or S4),
%         otherwise with first N2 or N3 or S4
%  N1_XR  starting with first N1 directly followed by either N2, N3, (S4), or R,
%         otherwise with first N2 or N3 or (S4) or R
%  NR     starting with first one of N2, N3 or (S4)
%  XR     starting with first one of N2, N3 or (S4) or R
%
%  For example N1_NR, sleep onset epoch in brackets:
%  [N1] N1 N2, but not [N1] N1 X N2, where X is not NR

% set the defaults
cfg.sleeponsetdef  = ft_getopt(cfg, 'sleeponsetdef', 'N1_XR');

hasLightsOff = false;
lightsOffMoment = 0;
if isfield(scoring, 'lightsoff')
    hasLightsOff = true;
    lightsOffMoment = scoring.lightsoff;
else
    ft_warning('The lights off moment was not provided in the scoring structure.\n The beginning of the scoring is thus assumed as lights off.');
end

epochs = cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0);
hypnStages = [cellfun(@sleepStage2str,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt2,epochs,'UniformOutput',0)];

onsetnumber = -1;

switch cfg.sleeponsetdef
    case 'NR'
        for iOnset = 1:numel(scoring.epochs)
            if strcmp(hypnStages(iOnset,3),'NR') && ((iOnset-1)*scoring.epochlength >= lightsOffMoment)
                onsetnumber = iOnset;
                break;
            end
        end
        
    case 'XR'
        for iOnset = 1:numel(scoring.epochs)
            if (strcmp(hypnStages(iOnset,3),'NR') ||  strcmp(hypnStages(iOnset,3),'R')) && ((iOnset-1)*scoring.epochlength >= lightsOffMoment)
                onsetnumber = iOnset;
                break;
            end
        end
        
    case {'N1' 'N1_NR' 'N1_XR'}
        
        consecN1 = 0;
        hasN1 = logical(0);
        for iOnset = 1:numel(scoring.epochs)
            if strcmp(hypnStages(iOnset,1),'N1') && ((iOnset-1)*scoring.epochlength >= lightsOffMoment)
                hasN1 = logical(1);
                consecN1 = consecN1 + 1;
                if ((onsetnumber+consecN1) ~= iOnset)
                    onsetnumber = iOnset;
                    consecN1 = 0;
                end
                if strcmp(cfg.sleeponsetdef,'N1')
                    break;
                end
            elseif ( strcmp(cfg.sleeponsetdef,'N1_XR') && (strcmp(hypnStages(iOnset,3),'NR') || strcmp(hypnStages(iOnset,3),'R')) ) ...
                    || ( strcmp(cfg.sleeponsetdef,'N1_NR') && (strcmp(hypnStages(iOnset,3),'NR')) ) ...
                    && ((iOnset-1)*scoring.epochlength >= lightsOffMoment)
                if ~hasN1
                    onsetnumber = iOnset;
                end
                break;
            else
                consecN1 = 0;
                hasN1 = logical(0);
            end
        end
end

preoffsetnumber = max(find(strcmp(hypnStages(:,1),'N1') | strcmp(hypnStages(:,3),'NR') | strcmp(hypnStages(:,3),'R') | strcmp(hypnStages(:,3),'MT')));

if isempty(preoffsetnumber);
    ft_warning('could not identify a sleep offset.\nAssume that sleep did not occur and thus not end.')
    preoffsetnumber = numel(scoring.epochs);
end

if onsetnumber == -1;
    ft_warning('could not identify a sleep onset with the current sleep onset definitions.\nAssume that sleep did not occur.')
    onsetnumber = numel(scoring.epochs);
    onsetepoch = '';
else
    onsetepoch = scoring.epochs{onsetnumber};
end

end

