function scoring = st_scoring_change_epochlength(cfg, scoring)

% ST_SCORING_CHANGE_EPOCH_LENGTH changes the epoch length (e.g 30 s) to a
% desired epoch length with interploation/extrapolation of epochs
%
% Use as
%   [scoring] = st_scoring_change_epochlength(cfg, scoring)
%
% Suggested configuration parameters are:
%   cfg.new_epochlength_seconds  = the desired epoch length 
%                                  (default is no change)
%
% See also ST_READ_SCORING, ST_SCORINGCONVERT

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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
st_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar scoring
ft_preamble provenance scoring
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set defaults
cfg.new_epochlength_seconds = ft_getopt(cfg, 'new_epochlength_seconds', scoring.epochlength);
                                    
fprintf([functionname ' function initialized\n']);

%scoring = scorings{1};
nEpochs = numel(scoring.epochs);
if scoring.epochlength == cfg.new_epochlength_seconds
    scoring_interpolated = scoring;
    return
end
interpolationfactor = scoring.epochlength/cfg.new_epochlength_seconds;

if interpolationfactor ~= round(interpolationfactor)
    ft_warning('the interpolation resulted in a non-integer epoch length!')
end

if nEpochs > 1
    res_cycle = [];
    res_cycle.ori = 'st_sleepcycles';
    res_cycle.type = 'cycle';
    res_cycle.cfg = [];
    res_cycle.table = table(1,1,nEpochs,nEpochs,'VariableNames',{'cycle','startepoch','endepoch','durationepochs'});
    cfg.newcycledurations = interpolationfactor*nEpochs;
    [scoring_interpolated] = st_scoringnormcycle(cfg, scoring, res_cycle);
    scoring_interpolated.epochlength = scoring_interpolated.epochlength * numel(scoring.epochs)/numel(scoring_interpolated.epochs);
    scoring = scoring_interpolated;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous scoring
ft_postamble provenance scoring
ft_postamble history    scoring
ft_postamble savevar    scoring


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)

end