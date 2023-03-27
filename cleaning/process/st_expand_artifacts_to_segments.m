function cfg_artifacts=st_expand_artifacts_to_segments(cfg_artifacts)

% ST_EXPAND_ARTIFACTS_TO_SEGMENTS expands artifacts from bad channels across segments
%
% Expand artifacts across all segments of channel when proportion of unrejected artifacual segments >= badchannelthresh
%
% Use as:
%     cfg=st_expand_artifacts_to_segments(cfg)
%
% Required configuration parameters:
%     cfg.grid      = structure containing segment-based artifact grids
%
% Optional configuration parameters (subfield grid):
%     cfg.badchannelthresh = minimum proportion of unrejected segments required to label all segments of a channel as artifact (default: 0.5)
%
% Output:
%     cfg = artifact configuration with added artifact grids:
%     - cfg.grid.artifact_grid_segment_expanded
%     - cfg.grid.artifact_grid_segment_expansion
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
numSeg=cfg_grid.segment_number;


%starting point is channel-expanded grid (= basic artifact + spatial expansion) 
artifact_grid_channel_expanded=cfg_grid.artifact_grid_channel_expanded;

%whatever has already been rejected needs to be excluded from consideration
rejection_grid=cfg_grid.reject_grid;


%take from cfg, otherwise default value
cfg_artifacts.badchannelthresh = ft_getopt(cfg_artifacts, 'badchannelthresh', 0.5);%propotion of bad segments to label channel for full interpolation

%segment-expanded grid:
%1) take the accepted/unrejected part of the artifact grid (using first chan of rejection grid)
%2) average across columns -> proportion of segments each channel is bad
%3) find channels above threshold, and repmat to expand to all segments
%(including rejected!)
expanded_seg_mat=artifact_grid_channel_expanded | repmat(mean(artifact_grid_channel_expanded(:,~rejection_grid(1,:)),2)>=cfg_artifacts.badchannelthresh,[1 numSeg]);

%expansion-only: disregard the original artifact grid
expansion_seg_mat=expanded_seg_mat & ~artifact_grid_channel_expanded;

%assign to cfg
cfg_grid.artifact_grid_segment_expanded=expanded_seg_mat;
cfg_grid.artifact_grid_segment_expansion=expansion_seg_mat;

cfg_artifacts.grid=cfg_grid;

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
