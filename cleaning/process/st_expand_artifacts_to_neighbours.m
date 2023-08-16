function cfg_artifacts=st_expand_artifacts_to_neighbours(cfg_artifacts)

% ST_EXPAND_ARTIFACTS_TO_NEIGHBOURS expands grid-based artifacts to neighboring channels
%
%   Non-artifactual channels are labeled artifactual when proportion of artifactual neighbouring channels >= channelexpandthresh.
%   Lower values are more aggressive (a few artifactual neighbours will turn a clean channel artifactual [0: expand even if NO neighbours are artifacts (= always set to artifact)]),
%   higher values are more conservative (many artifactual neighbours are needed to turn a clean channel artifactual [1: only expand if ALL neighbours are artifacts ; >1 (or Inf): never expand]).
%   Note: expansion behavior depends on both cfg.grid.channelexpandthresh and cfg.neighbours: a sensible neighbourhood structure is required.
%
% Use as:
%     cfg=st_expand_artifacts_to_neighbours(cfg)
%
% Required configuration parameters (automatically supplied by preceding function ST_EVENT_TABLES_TO_ARRAY):
%     cfg.grid      = structure containing segment-based artifact grids
%     cfg.neighbours = neighbourhod structure
%
% Optional configuration parameters:
%     cfg.channelexpandthresh = minimum proportion of artifactual neighbors required to expand artifacts (default: depending on number of channels, see below)
%      
% Output:
%     cfg = artifact configuration with added artifact grids:
%     - cfg.grid.artifact_grid_channel_expanded
%     - cfg.grid.artifact_grid_channel_expansion
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
ft_checkconfig(cfg_artifacts,'required',{'neighbours'});
ft_checkconfig(cfg_artifacts.artifacts,'required',{'grid'});

cfg_grid=cfg_artifacts.artifacts.grid;


artifact_grid=cfg_grid.artifact_grid_merged;

numChan=cfg_grid.channel_number;
numSeg=cfg_grid.segment_number;

%defaults for expanding artifacts to neighbors:
%become less aggressive with fewer channels
if numChan>128
    channelexpandthresh=0.25;
elseif numChan>64
    channelexpandthresh=0.5;
elseif numChan>10
    channelexpandthresh=0.75;
else %with very few channels, do not expand artifacts to neighbors
    channelexpandthresh=inf;
end

cfg_artifacts.channelexpandthresh=ft_getopt(cfg_artifacts,'channelexpandthresh',channelexpandthresh);


%get connectivity matrix (using neighbours)
cfg_artifacts.channel=cfg_artifacts.elec.label;
connectivity=channelconnectivity(cfg_artifacts);

%setup channel-expanded grid
artifact_grid_channel_expanded=artifact_grid;
for seg_i=1:numSeg

    %connected to artifact
    connected_and_artifact=double(connectivity & repmat(artifact_grid_channel_expanded(:,seg_i)',[numChan 1]));
    connected_and_artifact(~connectivity)=nan;
    artifact_grid_channel_expanded(:,seg_i)=artifact_grid_channel_expanded(:,seg_i) | nanmean(connected_and_artifact,2)>=cfg_artifacts.channelexpandthresh;
end

%also keep channel-expansion only
artifact_grid_channel_expansion=artifact_grid_channel_expanded & ~artifact_grid;

%assign to cfg
cfg_grid.artifact_grid_channel_expanded=artifact_grid_channel_expanded;
cfg_grid.artifact_grid_channel_expansion=artifact_grid_channel_expansion;

cfg_artifacts.artifacts.grid=cfg_grid;

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)