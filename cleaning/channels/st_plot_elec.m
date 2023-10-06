function fh=st_plot_elec(cfg)

% ST_PLOT_ELEC plots a blank topography with channel labels.
% 
% Use as:
%     st_plot_elec(cfg)
%     fh=st_plot_elec(cfg)
% 
% Required configuration parameters:
%     cfg.elec      = structure, containing channel information (FieldTrip format)
%
% Optional output:
%     fh        = handle to figure
%
% See also FT_READ_SENS

% Copyright (C) 2022-, Roy Cox, Frederik D. Weber
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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

%-----core function start---
%---input checks and defaults----
ft_checkconfig(cfg,'required','elec');
cfg.skipscale = ft_getopt(cfg,'skipscale', 'yes');
cfg.skipcomnt = ft_getopt(cfg,'skipcomnt','yes');

fprintf([functionname ' function initialized\n']);

%create layout from elec
layout = ft_prepare_layout(cfg);

%plot
fh=figure;
ft_plot_layout(layout,'label','yes','point','no','box','no'); %labels only
set(fh,'color','w');
%----core function end--


% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)