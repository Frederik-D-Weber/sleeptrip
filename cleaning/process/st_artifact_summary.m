function cfg_artifacts=st_artifact_summary(cfg_artifacts)

% ST_ARTIFACT_SUMMARY calculates a summary table based on information in
% cfg_artifacts.artifacts.grid. If a scoring is present, details are also provided by sleep
% stage.
%
% Use as:
%     cfg_artifacts=st_artifact_summary(cfg_artifacts)
%
% Required configuration parameters:
%     cfg_artifacts.artifacts.grid      = structure containing segment-based artifact grids
%
%Output:
%     cfg_artifacts = artifact configuration with added subfield:
%     - cfg_artifacts.artifacts.summary
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

hasScoring=false;
if isfield(cfg_artifacts,'scoring')
    scoring=cfg_artifacts.scoring; %may not exist
    hasScoring=true;
end

%repair and rejection summary
reject_grid=cfg_grid.reject_grid;
% repair_grid=cfg_grid.repair_grid;


% concatenate all artifact grids
artifact_labels=[cfg_grid.label 'BASIC' 'spatial expansion' 'temporal expansion' 'ANY'];
artifact_grids=cat(1,...
    cfg_grid.artifact_grid_by_type,...
    permute(cfg_grid.artifact_grid_merged,[3 1 2]),...
    permute(cfg_grid.artifact_grid_channel_expansion,[3 1 2]),...
    permute(cfg_grid.artifact_grid_segment_expansion,[3 1 2]),...
    permute(cfg_grid.artifact_grid_segment_expanded,[3 1 2])); %this is the final grid containing ALL artifacts plus ALL expansions

%same for other grids
global_labels={'REJECT','REPAIR'};
global_grids=cat(1,permute(cfg_grid.reject_grid,[3 1 2]),...
    permute(cfg_grid.repair_grid,[3 1 2]));

%pool
artifact_labels=[artifact_labels global_labels];
artifact_grids=cat(1,artifact_grids,global_grids);

%percentage_and_minutes= @(x) [100*mean(x,[2 3]) mean(x,[2 3])*size(x,3)*cfg_grid.segment_length/60];
calc_percentage= @(x) 100*mean(x,[ 2 3] );

%glob_pct=calc_percentage(global_grids)
art_pct=calc_percentage(artifact_grids);
art_nonrejected_pct=calc_percentage(artifact_grids(:,:,~reject_grid(1,:)));

summaryTable=array2table([art_pct art_nonrejected_pct],'RowNames',artifact_labels,'VariableNames',{'artifact_data_global_perc','artifact_nonrejected_global_perc'});



%%
%get segment-wise sleep scores (length dependent on scoring length)
if hasScoring

    %get artfact-level scoring
    cfg_tmp=[];
    cfg_tmp.scoring=scoring;
    cfg_tmp.cfg_artifacts=cfg_artifacts;

    scoring_artifact_level=st_scoring_artifact_level(cfg_tmp);

    %extract epoch list:
    epochList=scoring_artifact_level.epochs;
    scoring_labels=scoring_artifact_level.label;

    for label_i=1:length(scoring_labels)

        currLabel=scoring_labels{label_i};

        %calculate stage-wise percentages (all of stage, and non-rejected
        %of stage)
        art_pct=calc_percentage(artifact_grids(:,:,strcmp(epochList,currLabel)));
        art_nonrejected_pct=calc_percentage(artifact_grids(:,:,strcmp(epochList,currLabel) & ~reject_grid(1,:)));

        %art_nonrejected_pct=calc_percentage(artifact_grids(:,:,strcmp(scoring_artifact_level.epochs(~reject_grid(1,:)),currLabel)));

        varName=strjoin({'artifact' 'data', currLabel,'perc'},'_');
        varName_2=strjoin({'artifact' 'nonrejected', currLabel,'perc'},'_');

        summaryTable=[summaryTable array2table([art_pct art_nonrejected_pct],'RowNames',artifact_labels,'VariableNames',{varName varName_2})];

    end

    %scorings
    cfg_artifacts.scoring_artifact_level=scoring_artifact_level;

end


cfg_artifacts.artifacts.summary=summaryTable;

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
