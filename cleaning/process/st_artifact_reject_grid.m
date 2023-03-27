function cfg_artifacts=st_artifact_reject_grid(cfg_artifacts)

% ST_ARTIFACT_REJECT_GRID determines a rejection grid (i.e. segments to fully reject)
%
% Reject segment when proportion of artifactual chans >= segmentrejectthresh
% Lower values are more aggressive (a few artifactual channels will label a segment for rejection [0: always reject]),
% higher values are more conservative [1: only reject if all channels are artifactual, >1 (or Inf): never reject]
%
% Use as:
%     cfg=st_artifact_reject_grid(cfg)
%
% Required configuration parameters:
%     cfg.grid      = structure containing segment-based artifact grids
%
% Optional configuration parameters (subfield grid):
%     cfg.segmentrejectthresh = minimum proportion of artifactual channels required to label entire segment for rejection (default: depending on number of channels, see below)
%
% Output:
%     cfg = artifact configuration with added artifact grids:
%     - cfg.grid.reject_grid
%     - cfg.grid.accept_grid
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
ft_checkconfig(cfg_artifacts,'required',{'grid'});
cfg_grid=cfg_artifacts.grid;

%use the channel-expanded grid (= basic artifact + spatial expansion)
artifact_grid=cfg_grid.artifact_grid_channel_expanded;

numChan=size(artifact_grid,1);

%segment rejection:
%become less aggressive with fewer channels
if numChan>128
    segmentrejectthresh=0.25; %proportion of artifactual channels to label segment for rejection
elseif numChan>10
    segmentrejectthresh=0.5;
else
    segmentrejectthresh=0.75;
end

%since we can't repair with too few channels, set threshold such that rejection occurs if any channel is bad
if numChan<3 
    segmentrejectthresh=1/numChan;
end

%take from cfg, otherwise default set above
cfg_artifacts.segmentrejectthresh = ft_getopt(cfg_artifacts, 'segmentrejectthresh', segmentrejectthresh);

reject_grid=repmat(mean(artifact_grid,1)>=cfg_artifacts.segmentrejectthresh,[numChan 1]);
accept_grid=~reject_grid;

%assign to cfg
cfg_grid.reject_grid=reject_grid;
cfg_grid.accept_grid=accept_grid;

cfg_artifacts.grid=cfg_grid;

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)