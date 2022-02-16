function [onsetnumber lastscoredsleepstagenumber onsetepoch lastscoredsleepstage allowedsleeponsetbeforesleepopon allowedsleepaftersleepopoff] = st_sleeponset(cfg,scoring)
% 
% ST_SLEEPONSET determines the sleep onset of a sleep scoring
% Use as
%   [onsetnumber, lastscoredsleepstagenumber, onsetepoch, lastscoredsleepstage, allowedsleeponsetbeforesleepopon allowedsleepaftersleepopoff] = st_sleeponset(cfg,scoring)
%   [onsetnumber, lastscoredsleepstagenumber, onsetepoch, lastscoredsleepstage, allowedsleeponsetbeforesleepopon] = st_sleeponset(cfg,scoring)
%   [onsetnumber, lastscoredsleepstagenumber, onsetepoch, lastscoredsleepstage] = st_sleeponset(cfg,scoring)
%   [onsetnumber, lastscoredsleepstagenumber] = st_sleeponset(cfg,scoring)
%   [onsetnumber] = st_sleeponset(cfg,scoring)
%
% Configutation parameter can be empty, e.g. cfg = []
%
% Optional configuration parameters are
%   cfg.sleeponsetdef  = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                        'NR' or 'N2R' or 'XR' or 'AASM' or 'X2R' or 
%                        'N2' or 'N3' or 'SWS' or 'S4' or 'R', 
%                        see below for details (default = 'N1_XR')
%   cfg.allowsleeponsetbeforesleepopon = srting, if possible, allow sleep onset before sleep
%                        opportunity (or lights off moment if former is not present) 
%                        either 'yes' or 'no' (default = 'no')
%   cfg.allowsleepaftersleepopoff = srting, if possible, allow sleep (offset, i.e. end of sleep) after sleep
%                        opportunity (or lights on moment if former is not present) 
%                        either 'yes' or 'no' see ST_SLEEPONSET for details (default = 'yes')
%
%  Here are the possible sleep onset definitions, where NR referes to
%  non-REM, R to REM and XR to any non-REM or REM sleep stage.
%  N1     starting with first N1 is sleep onset, nothing else
%  N2     starting with first N2 is sleep onset, nothing else
%  N3     starting with first N3 is sleep onset, nothing else
%  S4     starting with first S4 is sleep onset, nothing else
%  SWS    starting with first SWS is sleep onset, i.e. S3/N3 or S4, nothing else
%  R      starting with first R is sleep onset, nothing else
%  N1_NR  starting with first N1 followed directly by either N2, N3 (or S4),
%         otherwise with first N2 or N3 or S4
%  N1_XR  starting with first N1 directly followed by either N2, N3, (S4), or R,
%         otherwise with first N2 or N3 or (S4) or R
%  NR      starting with first one of N1, N2, N3 or (S4)
%  N2R     starting with first one of N2, N3 or (S4)
%  XR     starting with first one of N1, N2, N3 or (S4) or R
%  AASM   starting with first one of N1, N2, N3 or (S4) or R
%  X2R     starting with first one of N2, N3 or (S4) or R
%
%  For example N1_NR, sleep onset epoch in brackets:
%  [N1] N1 N2, but not [N1] N1 X N2, where X is not NR

% set the defaults
cfg.sleeponsetdef  = upper(ft_getopt(cfg, 'sleeponsetdef', 'N1_XR'));
cfg.allowsleeponsetbeforesleepopon  = ft_getopt(cfg, 'allowsleeponsetbeforesleepopon', 'no');
cfg.allowsleepaftersleepopoff  = ft_getopt(cfg, 'allowsleepaftersleepopoff', 'no');




hasLightsOff = false;
lightsOffMoment = 0;
if isfield(scoring, 'lightsoff')
    hasLightsOff = true;
    lightsOffMoment = scoring.lightsoff;
else
    ft_warning('The lights off moment was not provided in the scoring structure.\n The beginning of the scoring is thus assumed as lights off.');
end


hasLightsOn = false;
lightsOnMoment = NaN;
if isfield(scoring, 'lightson')
    hasLightsOn = true;
    lightsOnMoment = scoring.lightson;
else
    ft_warning('The lights on moment was not provided in the scoring structure.\n The last sleep stage is assumed to best match this.');
end


hasSleepOpportunityOn = false;
sleepOpportunityOnMoment = 0;
if isfield(scoring, 'sleepopon') 
    if ~isnan(scoring.sleepopon)
        hasSleepOpportunityOn = true;
        sleepOpportunityOnMoment = scoring.sleepopon;
    else
        ft_warning('The sleep opportunity onset moment was NaN in the scoring structure.\n The beginning of the scoring is thus assumed as sleep opportunity onset, but sleep onset will be NaN.');
    end
else
    if hasLightsOff && ~isnan(lightsOffMoment)
        sleepOpportunityOnMoment = lightsOffMoment;
        hasSleepOpportunityOn = true;
    	ft_warning('The sleep opportunity onset moment was not provided in the scoring structure.\n The lights off moment is used instead');
    else
    	ft_warning('The sleep opportunity onset moment was not provided in the scoring structure.\n The beginning of the scoring is thus assumed as sleep opportunity onset.');
    end
end

hasSleepOpportunityOff = false;
sleepOpportunityOffMoment = NaN;
if isfield(scoring, 'sleepopoff')
    if ~isnan(scoring.sleepopon)
        hasSleepOpportunityOff = true;
        sleepOpportunityOffMoment = scoring.sleepopoff;
    else
        ft_warning('The sleep opportunity offset moment was NaN in the scoring structure.\n The end of the scoring is thus assumed as sleep opportunity onset, but sleep opportunity window and period related measures will be NaN.');
    end
else
    if hasLightsOn && ~isnan(lightsOnMoment)
        sleepOpportunityOffMoment = lightsOnMoment;
        hasSleepOpportunityOff = true;
    	ft_warning('The sleep opportunity off moment was not provided in the scoring structure.\n The lights on moment is used instead');
    else
    	ft_warning('The sleep opportunity off moment was not provided in the scoring structure.\n The last sleep stage is assumed to best match this.');
    end
end


epochs = cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0);
hypnStages = [cellfun(@sleepStage2str,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt2,epochs,'UniformOutput',0)...
    cellfun(@sleepStage2str_alt3,epochs,'UniformOutput',0)];

onsetnumber = -1;


allowedsleeponsetbeforesleepopon = false;
if hasSleepOpportunityOn
        for iOnset = 1:numel(scoring.epochs)
            if strcmp(hypnStages(iOnset,4),'S') 
                break;
            end
        end  
        if (iOnset-1)*scoring.epochlength < sleepOpportunityOnMoment
                if istrue(cfg.allowsleeponsetbeforesleepopon)
                    allowedsleeponsetbeforesleepopon = true;
                    ft_warning('There were sleep stages scored at epoch %d BEFORE the sleep opportunity onset (which might have defaulted to ligths off moment) at %f s!\n BUT sleep onset is allowed to start before sleep opportunity onset!',iOnset,sleepOpportunityOnMoment);
                else
                	ft_warning('There were sleep stages scored at epoch %d BEFORE the sleep opportunity onset (which might have defaulted to ligths off moment) at %f s!\n The epoch after sleep opportunity onset is thus assumed as sleep onset. or use the cfg.allowsleeponsetbeforesleepopon = ''yes'' to ignore this case!',iOnset,sleepOpportunityOnMoment);
                end
        end
end


switch cfg.sleeponsetdef
     case 'NR'
        for iOnset = 1:numel(scoring.epochs)
            if (strcmp(hypnStages(iOnset,1),'N1') || strcmp(hypnStages(iOnset,3),'NR')) && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end  
     case 'N2'
        for iOnset = 1:numel(scoring.epochs)
            if ( strcmp(hypnStages(iOnset,1),'N2') ) && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end 
      case 'N3'
        for iOnset = 1:numel(scoring.epochs)
            if ( strcmp(hypnStages(iOnset,1),'N3') ) && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end  
     case 'S4'
        for iOnset = 1:numel(scoring.epochs)
            if ( strcmp(hypnStages(iOnset,1),'S4') ) && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end 
    case 'R'
        for iOnset = 1:numel(scoring.epochs)
            if ( strcmp(hypnStages(iOnset,1),'R') ) && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end  
    case 'SWS'
        for iOnset = 1:numel(scoring.epochs)
            if ( strcmp(hypnStages(iOnset,1),'N3') || strcmp(hypnStages(iOnset,1),'S4') ) && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end  
    case 'N2R'
        for iOnset = 1:numel(scoring.epochs)
            if strcmp(hypnStages(iOnset,3),'NR') && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end  

     case 'X2R'
        for iOnset = 1:numel(scoring.epochs)
            if (strcmp(hypnStages(iOnset,3),'NR') ||  strcmp(hypnStages(iOnset,3),'R')) && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end      
        
    case {'XR' 'AASM'}
        for iOnset = 1:numel(scoring.epochs)
            if (strcmp(hypnStages(iOnset,1),'N1') || strcmp(hypnStages(iOnset,3),'NR') ||  strcmp(hypnStages(iOnset,3),'R')) && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment)  || allowedsleeponsetbeforesleepopon)
                onsetnumber = iOnset;
                break;
            end
        end
        
    case {'N1' 'N1_NR' 'N1_XR'}
        
        consecN1 = 0;
        hasN1 = logical(0);
        for iOnset = 1:numel(scoring.epochs)
            if strcmp(hypnStages(iOnset,1),'N1') && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
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
                    && (((iOnset-1)*scoring.epochlength >= sleepOpportunityOnMoment) || allowedsleeponsetbeforesleepopon)
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



lastscoredsleepstagenumber = max(find(strcmp(hypnStages(:,1),'N1') | strcmp(hypnStages(:,3),'NR') | strcmp(hypnStages(:,3),'R') | strcmp(hypnStages(:,3),'MT')));

allowedsleepaftersleepopoff = false;
if hasSleepOpportunityOff
%         if (scoring.epochlength*numel(scoring.epochs)) > sleepOpportunityOffMoment
%             
%         else
        if lastscoredsleepstagenumber*scoring.epochlength > sleepOpportunityOffMoment
                if istrue(cfg.allowsleepaftersleepopoff)
                    allowedsleepaftersleepopoff = true;
                    ft_warning('There were sleep stages scored at epoch %d AFTER the sleep opportunity off (which might have defaulted to ligths on moment) at %f s!\n BUT sleep off is allowed to end after sleep opportunity off!',lastscoredsleepstagenumber,sleepOpportunityOffMoment);
                else
                    temp_sleepopoff_epoch = min(numel(hypnStages(:,1)),ceil(sleepOpportunityOffMoment/scoring.epochlength));
                    lastscoredsleepstagenumber_sleep_opoff = max(find(strcmp(hypnStages(1:temp_sleepopoff_epoch,1),'N1') | strcmp(hypnStages(1:temp_sleepopoff_epoch,3),'NR') | strcmp(hypnStages(1:temp_sleepopoff_epoch,3),'R') | strcmp(hypnStages(1:temp_sleepopoff_epoch,3),'MT')));
                    lastscoredsleepstagenumber = min(lastscoredsleepstagenumber,lastscoredsleepstagenumber_sleep_opoff);
                	ft_warning('There were sleep stages scored at epoch %d AFTER the sleep opportunity off (which might have defaulted to ligths on moment) at %f s!\n The epoch in which the sleep opportunity off ends is thus assumed as sleep offset. or use the cfg.allowsleepaftersleepopoff = ''yes'' to ignore this case!',iOnset,sleepOpportunityOnMoment);
                end
        end
%         end
        

else
    
end


if isempty(lastscoredsleepstagenumber) || isnan(lastscoredsleepstagenumber)
    ft_warning('could not identify a sleep offset.\nAssume that sleep did not occur and thus not end.')
    lastscoredsleepstagenumber = numel(scoring.epochs);
end

if (lastscoredsleepstagenumber > 0) && (lastscoredsleepstagenumber <= numel(scoring.epochs))
lastscoredsleepstage = scoring.epochs{lastscoredsleepstagenumber};
else
  lastscoredsleepstage = '';  
end

if onsetnumber == -1;
    ft_warning('could not identify a sleep onset with the current sleep onset definitions.\nAssume that sleep did not occur.')
    onsetnumber = numel(scoring.epochs);
    onsetepoch = '';
else
    onsetepoch = scoring.epochs{onsetnumber};
end

end

