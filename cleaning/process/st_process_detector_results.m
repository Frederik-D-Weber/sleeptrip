function cfg_artifacts=st_process_detector_results(cfg_artifacts)

% ST_PROCESS_DETECTOR_RESULTS processes artifacts (as detected by ST_RUN_DETECTOR_SET) as follows:
% - convert continuous artifacts (event tables) to segment-based artifact grids (ST_EVENT_TABLES_TO_ARRAY)
% - expand artifacts to neighboring channels (ST_EXPAND_ARTIFACTS_TO_NEIGHBOURS)
% - determine rejection grid (ST_ARTIFACT_REJECT_GRID)
% - determine bad channels (expanding artifacts across all segments) (ST_EXPAND_ARTIFACTS_TO_SEGMENTS)
% - determine repair matrix (ST_REPAIR_GRID)
%- calculate summary table (ST_ARTIFACT_SUMMARY)
% - (optionally) turn all grids into event tables (for export) (ST_ARTIFACT_SEGMENT_TO_CONTINUOUS)
%
% Use as:
%     cfg=st_process_detector_results(cfg)
%
% Required configuration parameters (all automatically supplied by ST_RUN_DETECTOR_SET):
%     cfg.continuous      = structure containing continuous artifact event tables
%     cfg.detector_set  = cfg containing details of individual detectors
%     cfg.data = data structure
%     cfg.elec = elec structure
%
% Optional configuration parameters:
%     [main]
%     cfg.neighbours = neighbourhood structure to be used by ST_EXPAND_ARTIFACTS_TO_NEIGHBOURS (default taken from ST_GET_DEFAULT_NEIGHBOURS)
%     cfg.keepeventtables = string (default: 'no')
%     cfg.merge_detectors = cell array of strings, with names of each detector to include for computation of artifact grid. (default: 'all' [string])
%     cfg.segment_length = length of grid segments in seconds (default: 5)
%     cfg.channelexpandthresh = minimum proportion of artifactual neighbors required to expand artifacts (default: depending on number of channels, see ST_EXPAND_ARTIFACTS_TO_NEIGHBOURS)
%     cfg.segmentrejectthresh = minimum proportion of channels required to label segment as rejected (default: depending on number of channels, see ST_ARTIFACT_REJECT_GRID)
%     cfg.badchannelthresh = minimum proportion of unrejected segments required to label all segments of a channel as artifact (default: 0.5, see ST_EXPAND_ARTIFACTS_TO_SEGMENTS)
%
% Output:
%     cfg = configuration containing channelwise artifact information
%
% See also ST_GET_DEFAULT_DETECTOR_SET ST_RUN_DETECTOR_SET

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
%ft_preamble init
% ft_preamble debug
% ft_preamble loadvar data
% ft_preamble provenance data
% ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
% if ft_abort
%   return
% end

%-----core function start---
%---input checks and defaults----
ft_checkconfig(cfg_artifacts,'required',{'continuous','detector_set','data','elec'});
cfg_artifacts.keepeventtables  = ft_getopt(cfg_artifacts, 'keepeventtables', 'yes');

fprintf([functionname ' function initialized\n'])


if ~isfield(cfg_artifacts,'neighbours')
    %get default neighbours from elec
    tmp_cfg=[];
    tmp_cfg.elec=cfg_artifacts.elec;
    neighbours=st_get_default_neighbours(tmp_cfg); %calculate default neighbours

    cfg_artifacts.neighbours= ft_getopt(cfg_artifacts,'neighbours',neighbours);

end



%%

%--------PROCESS EVENTS-------

%--------------1: convert continuous artifacts (requested types only) to basic artifact_grid (ch x seg)

%get artifact grids (both merged and separate by detector type)
cfg_artifacts=st_event_tables_to_array(cfg_artifacts);

%-------------2: channel expanded artifacts (based on neighbors)
%channel-expand artifacts based on neighborhood structure
cfg_artifacts=st_expand_artifacts_to_neighbours(cfg_artifacts);

%-------3: segments to fully reject
cfg_artifacts=st_artifact_reject_grid(cfg_artifacts);

%-----------4: channels to fully repair (ignoring rejected segments)
cfg_artifacts=st_expand_artifacts_to_segments(cfg_artifacts);

%--------------5: repair matrix (= original or time-interpolated artifacts, minus
cfg_artifacts=st_repair_grid(cfg_artifacts);



%-- 6 calculate summary table
cfg_artifacts=st_artifact_summary(cfg_artifacts);


if strcmp(cfg_artifacts.keepeventtables,'yes')
    %turn all grids into to continuous event tables
    cfg_artifacts=st_artifact_segment_to_continuous(cfg_artifacts);

    %remove the continuous event tables (memory intensive)
elseif strcmp(cfg_artifacts.keepeventtables,'no')
    cfg_artifacts=removefields(cfg_artifacts,{'continuous','continuous_from_grid'});
end

% do the general cleanup and bookkeeping at the end of the function
% ft_postamble debug
% ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)