function cfg_artifacts=st_repair_grid(cfg_artifacts)

% ST_REPAIR_GRID determines the repair grid/matrix
%
% Use as:
%     cfg=st_repair_grid(cfg)
%
% Required configuration parameters:
%     cfg.grid      = structure containing segment-based artifact grids
%
% Output:
%     cfg = artifact configuration with added grids:
%     - cfg.grid.repair_grid
%
% See also ST_PROCESS_DETECTOR_RESULTS ST_GET_DEFAULT_DETECTOR_SET ST_RUN_DETECTOR_SET

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

%---input checks and defaults----
ft_checkconfig(cfg_artifacts.artifacts,'required',{'grid'});
cfg_grid=cfg_artifacts.artifacts.grid;

%determine repair matrix (= original or channel-expanded or segment-expanded artifacts, minus rejected)
repair_grid=(cfg_grid.artifact_grid_merged | cfg_grid.artifact_grid_channel_expansion | cfg_grid.artifact_grid_segment_expansion) & ~cfg_grid.reject_grid;

numChan=cfg_grid.channel_number;
if numChan<3 %meaningless to repair
    repair_grid(:)=false;
end

cfg_grid.repair_grid=repair_grid;

cfg_artifacts.artifacts.grid=cfg_grid;

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)