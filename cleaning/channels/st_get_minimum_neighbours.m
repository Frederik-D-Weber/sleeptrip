function neighbours=st_get_minimum_neighbours(cfg)

% ST_GET_MINIMUM_NEIGHBOURS creates a neighborhood structure such that each channel has at least a specified number of neighbors.
%
% Use as:
%     neighbours=st_get_minimum_neighbors(cfg)
%
% Required configuration parameters:
%     cfg.elec      = FieldTrip elec/sens structure, containing channel information
%     cfg.minimumneighbours = minimum number of neighbors each channel should have
%
% Output:
%     neighbours        = FieldTrip neighbourhood structure
%
% See also ST_PLOT_NEIGHBOURS FT_PREPARE_NEIGHBOURS

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
ft_checkconfig(cfg,'required',{'elec','minimumneighbours'});
fprintf([functionname ' function initialized\n']);

%find neighbor distance yielding a minimum number of neighbors
minNumNeighb=cfg.minimumneighbours;
elec=cfg.elec;

chanpos=elec.chanpos;
numChan=length(elec.label);

dist = zeros(numChan,numChan);
for i=1:numChan
    dist(i,:) = sqrt(sum((chanpos(1:numChan,:) - repmat(chanpos(i,:), numChan, 1)).^2,2))';
end
sortDist=sort(dist,1);
sortDist(1,:)=[];%remove first row of zeros
maxDistToNeighb=max(sortDist,[],2);
if minNumNeighb<=numChan-2
    neighbDist=maxDistToNeighb(minNumNeighb+1);
else
    neighbDist=inf;
end

%set up neighborhood structure
cfg_nb=[];
cfg_nb.elec=elec;
cfg_nb.method = 'distance';
cfg_nb.neighbourdist=neighbDist;

neighbours=ft_prepare_neighbours(cfg_nb);

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