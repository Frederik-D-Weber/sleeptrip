function [scoring] = st_scoringconvert(cfg, scoring)

% ST_SCORINGCONVERT converts between scoring standards and hypnograms for 
% a structure provided by ST_READ_SCORING 
%
% Use as
%   [scoring] = st_scoringconvert(cfg, scoring)
%
% Required configuration parameters are:
%   cfg.to    = string, scoring standard to convert to either 
%               'aasm' for AASM
%               'rk' for Rechtschaffen&Kales or 
%               'number' or 'numbers' for conversion to {0, 1,  2,  3,  4,  5, -1} for {W, N1/S1,  N2/S2,  N3/S3,  N3/S4,  R, ?}
%               'custom' for providing a scoremap
% 
% in case cfg.to = 'custom' the following parameter also needs to be present
%   cfg.scoremap  = a structure must be provided (see below) 
% 
% in case cfg.to IS NOT 'custom' but the scoring.standard = 'custom' then
% either cfg.scoremap or if not also scoring.scoremap needs to be defined
%   then the optional parameter 
%    cfg.forceundefinedto = string, as a sleep stage to force the sleep stages to be converted into that are not part of the cfg.to 
%                           non-custom scoring (e.g. 'aasm' or 'rk' or 'numbers') 
%                           e.g. '?' or 'U'
%                           
%
%   Example scoremap sturcture, all three fields are required
%   scoremap = [];
%   scoremap.labelold  = {'0', '1',  '2',  '3',  '4',  '5', '8'};
%   scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', 'W'};
%   scoremap.unknown   = '?';
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

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

% set defaults


if strcmp(cfg.to,'custom') && ~isfield(cfg,'scoremap')
	ft_error('if the cfg.to paramter is set to ''custom'' it requires also a cfg.scoremap as parameter in the configuration.');
end

fprintf([functionname ' function initialized\n']);


iscustom = false;
scoremap = [];
scoremap.unknown   = '?';
switch scoring.standard
    case 'aasm'
        scoremap.labelold  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?', 'W'};
    case 'rk'
        scoremap.labelold  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?', 'MT'};
    case {'number', 'numbers'}
        scoremap.labelold  = {'0', '1', '2', '3', '4', '5', '-1', '8'};
    case 'custom'
        if ~isfield(cfg,'scoremap')
            if ~isfield(scoring,'scoremap')
                ft_error('the scoring was a ''custom'' label and thus it requires also a cfg.scoremap as parameter in the configuration.');
            else
                ft_warning('the scoring was a ''custom'' label and and the scoremap was taken from scoring.scoremap');
                cfg.scoremap = scoring.scoremap;
            end
        end
        scoremap.labelold  = cfg.scoremap.labelold;
        iscustom = true;
end


tononcustom = false;

switch cfg.to
    case 'aasm'
        scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?', 'W'};
        scoring.standard = 'aasm';
        tononcustom = true;
    case 'rk'
        scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?', 'MT'};
        scoring.standard = 'rk';
        tononcustom = true;
    case {'number', 'numbers'}
        scoremap.labelnew  = {'0', '1', '2', '3', '4', '5', '-1', '8'};
        scoring.standard = 'number';
        scoremap.unknown   = '-1';
        tononcustom = true;
    case 'custom'
    	if ~isfield(cfg,'scoremap')
            if ~isfield(scoring,'scoremap')
                ft_error('the scoring will be converted to a ''custom'' scoring and thus requires also a cfg.scoremap as parameter in the configuration.');
            else
                ft_warning('the scoring was a ''custom'' label and and the scoremap was taken from scoring.scoremap');
                cfg.scoremap = scoring.scoremap;
            end
        end
        scoremap.labelnew  = cfg.scoremap.labelnew;
        scoremap.unknown = cfg.scoremap.unknown;
        scoring.standard = 'custom';
    otherwise 
        ft_error(['the cfg.to option ' cfg.to ' is not known, please see the help of this function.']);
end

forceundefinedto = '';
if isfield(cfg,'forceundefinedto')
    forceundefinedto = cfg.forceundefinedto;
end

if iscustom && tononcustom
    iold_scoremap = ismember(cfg.scoremap.labelnew,scoremap.labelnew);
    if any(~iold_scoremap) && isempty(forceundefinedto)
        ft_error('The following scoring labels ''%s'' in the custom cfg.scoremap are not part of the standard in cfg.to = ''%s'' defined by the labels %s . Maybe you want to set cfg.forceundefinedto option to force a sleep stage (BUT THAT IS NOT RECOMMENDED).',strjoin(cfg.scoremap.labelnew(~iold_scoremap)),cfg.to,strjoin(unique(scoremap.labelnew)))
    end
    undefined_stages = {};
    if any(~iold_scoremap) && ~isempty(forceundefinedto)
        undefined_stages = cfg.scoremap.labelnew(~iold_scoremap);
        [iold inew] = ismember(scoring.epochs, undefined_stages);
        scoring.epochs(iold) = {forceundefinedto};
    end
    label = scoremap.labelnew';
    [iold inew] = ismember(scoring.epochs, cat(2,scoremap.labelnew, {forceundefinedto}));
    scoring.epochs(~iold) = {scoremap.unknown};
else
    epochs = repmat({scoremap.unknown},size(scoring.epochs));
    label = scoring.label;
    [labelold ia iItems] = unique(scoremap.labelold,'stable');
    usedItems = [];
    
    match_cum = logical(zeros(numel(scoring.epochs),1));
  
    
    for iLabel = 1:numel(scoremap.labelold)
        if ismember(iItems(iLabel),usedItems)
            continue
        end
        old = scoremap.labelold{iLabel};
        new = scoremap.labelnew{iLabel};
        match = cellfun(@(x) strcmp(x, old), scoring.epochs, 'UniformOutput', 1);
        epochs(match) = {new};
        match_cum = match_cum | logical(match(:));
        match = cellfun(@(x) strcmp(x, old), scoring.label, 'UniformOutput', 1);
        label(match) = {new};
        usedItems = [usedItems iItems(iLabel)];
    end
    
    if any(~match_cum)
        ft_warning('The sleep stages''%s'' in the original scoring were not covered in the scoremap and have thus been set to ''%s''',strjoin(unique(scoring.epochs(~match_cum))),scoremap.unknown)
    end
    
    scoring.epochs = epochs;
end
[scoring.label,ia,ic] = unique(label);

if isfield(scoring,'prob')
    prob = nan(numel(scoring.label),size(scoring.prob,2));
    for iLabel = 1:numel(scoring.label)
        lab = scoring.label{iLabel};
        prob(iLabel,:) = nansum(scoring.prob(ismember(label,lab),:),1);
    end
    scoring.prob = prob;
end




fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
