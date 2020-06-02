function [scoring] = st_scoringconvert(cfg, scoring)

% ST_SCORINGCONVERT converts between scoring standards and hypnograms for 
% a structure provided by ST_READ_SCORING 
%
% Use as
%   [cfg] = st_scoringconvert(cfg, scoring)
%
% Required configuration parameters are:
%   cfg.to    = string, scoring standard to convert to either 
%               'aasm' for AASM
%               'rk' for Rechtschaffen&Kales or 
%               'numbers' for conversion to {0, 1,  2,  3,  4,  5, -1} for {W, N1/S1,  N2/S2,  N3/S3,  N3/S4,  R, ?}
%               'custom' for providing a scoremap
% 
% in case cfg.to = 'custom' the following parameter also needs to be present
%   cfg.scoremap  = a structure must be provided (see below) 
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

tic
memtic
st = dbstack;
functionname = st.name;
fprintf([functionname ' function started\n']);

% set defaults


if strcmp(cfg.to,'custom') && ~isfield(cfg,'scoremap')
	ft_error('if the cfg.to paramter is set to ''custom'' it requires also a cfg.scoremap as parameter in the configuration.');
end

fprintf([functionname ' function initialized\n']);


scoremap = [];
scoremap.unknown   = '?';
switch scoring.standard
    case 'aasm'
        scoremap.labelold  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?', 'W'};
    case 'rk'
        scoremap.labelold  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?', 'MT'};
    case 'number'
        scoremap.labelold  = {'0', '1', '2', '3', '4', '5', '-1', '8'};
    case 'custom'
        if ~isfield(cfg,'scoremap')
            ft_error('the scoring was a ''custom'' label and thus it requires also a cfg.scoremap as parameter in the configuration.');
        end
        scoremap.labelold  = cfg.scoremap.labelold;
end


switch cfg.to
    case 'aasm'
        scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?', 'W'};
        scoring.standard = 'aasm';
    case 'rk'
        scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?', 'MT'};
        scoring.standard = 'rk';
    case 'number'
        scoremap.labelnew  = {'0', '1', '2', '3', '4', '5', '-1', '8'};
        scoring.standard = 'number';
        scoremap.unknown   = '-1';
    case 'custom'
        scoremap.labelnew  = cfg.scoremap.labelnew;
        scoring.standard = 'custom';
    otherwise 
        ft_error(['the cfg.to option ' cfg.to ' is not known, please see the help of this function.']);
end



epochs = repmat({scoremap.unknown},size(scoring.epochs));
label = scoring.label;
for iLabel = 1:numel(scoremap.labelold)
    old = scoremap.labelold{iLabel};
    new = scoremap.labelnew{iLabel};
    match = cellfun(@(x) strcmp(x, old), scoring.epochs, 'UniformOutput', 1);
    epochs(match) = {new};
    match = cellfun(@(x) strcmp(x, old), scoring.label, 'UniformOutput', 1);
    label(match) = {new};
end
scoring.epochs = epochs;
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
toc
memtoc
end
