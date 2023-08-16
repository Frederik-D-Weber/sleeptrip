function cfg_artifacts=st_artifact_summary(cfg_artifacts)

% ST_ARTIFACT_SUMMARY calculates a summary table based on information in
% cfg.grid. If a scoring is present, details are also provided by sleep
% stage.
%
% Use as:
%     cfg=st_artifact_summary(cfg)
%
% Required configuration parameters:
%     cfg.grid      = structure containing segment-based artifact grids
%
%Output:
%     cfg = artifact configuration with added field:
%     - cfg.artifact_summary
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


% %rejection grid (global)
% reject_total_perc=100*mean(reject_grid(1,:)); %percentage of total data
% reject_total_min=sum(reject_grid(1,:))*cfg_grid.segment_length/60; %minutes of total data
%
% %repair grid as perc of unrejected data
% repair_unrejected_total_perc=100*mean(repair_grid(:,~reject_grid(1,:)),'all');
%
%
% %repair grid (excluding rejected) as perc of total data
% repair_data_total_perc=100*mean(repair_grid,'all'); %percentage of total data
%
% summaryTable=table(reject_total_perc,reject_total_min,repair_unrejected_total_perc,repair_data_total_perc);
%%
%get segment-wise sleep scores (length dependent on scoring length)
if hasScoring

    segment_length=cfg_grid.segment_length;
    numArtifactSegment=cfg_grid.segment_number;

    %rescore to match artifact segments
    tmp_cfg=[];
    tmp_cfg.new_epochlength_seconds=segment_length;
    scoring_artifact_level = st_scoring_change_epochlength(tmp_cfg, scoring);

    numScoringSegment=length(scoring_artifact_level.epochs);

    %artifact-level scoring may be longer or shorter than artifact segments:
    if numArtifactSegment<numScoringSegment % trim scoring if needed
        scoring_artifact_level.epochs=scoring_artifact_level.epochs(1:numArtifactSegment);
        scoring_artifact_level.excluded=scoring_artifact_level.excluded(1:numArtifactSegment);
        scoring_artifact_level.prob=scoring_artifact_level.prob(:,1:numArtifactSegment);
        scoring_artifact_level.numbers=scoring_artifact_level.numbers(1:numArtifactSegment);
    elseif numArtifactSegment>numScoringSegment %extend scoring by repeating last epoch info
        scoring_artifact_level.epochs(end+1:numArtifactSegment)=scoring_artifact_level.epochs(end);
        scoring_artifact_level.excluded(end+1:numArtifactSegment)= scoring_artifact_level.excluded(end);
        scoring_artifact_level.prob(:,end+1:numArtifactSegment)=repmat(scoring_artifact_level.prob(:,end),[1 numArtifactSegment-numScoringSegment]);
        scoring_artifact_level.numbers(end+1:numArtifactSegment)= scoring_artifact_level.numbers(end);
    end

    %extract epoch lists:
    %regular
    epochList=scoring_artifact_level.epochs;
%     %non_rejected only
%     epochList_nonrejected=epochList(~reject_grid(1,:));

    scoring_labels=scoring_artifact_level.label;

    for label_i=1:length(scoring_labels)

        currLabel=scoring_labels{label_i};

        %         %----rejection
        %         %percentage of stage
        %         reject_stage_perc=100*mean(reject_grid(1,strcmp(scoring_artifact_level_exclude.epochs,currLabel)));
        %         varName=strjoin({'reject' currLabel,'perc'},'_');
        %         summaryTable=[summaryTable table(reject_stage_perc,'VariableNames',{varName})];
        %
        %         %minutes of stage
        %         reject_stage_min=sum(reject_grid(1,strcmp(scoring_artifact_level_exclude.epochs,currLabel)))*segment_length/60;
        %         varName=strjoin({'reject' currLabel,'min'},'_');
        %         summaryTable=[summaryTable table(reject_stage_min,'VariableNames',{varName})];
        %
        %         %--repair (of non-rejected)
        %         repair_unrejected_stage_perc=100*mean(repair_grid(:,strcmp(scoring_artifact_level_exclude.epochs,currLabel) & ~scoring_artifact_level_exclude.excluded),'all');
        %         varName=strjoin({'repair' 'unrejected', currLabel,'perc'},'_');
        %         summaryTable=[summaryTable table(repair_unrejected_stage_perc,'VariableNames',{varName})];
        %
        %         %--repair (of total)
        %         repair_data_stage_perc=100*mean(repair_grid(:,strcmp(scoring_artifact_level_exclude.epochs,currLabel)),'all');
        %         varName=strjoin({'repair' 'data', currLabel,'perc'},'_');
        %         summaryTable=[summaryTable table(repair_unrejected_stage_perc,'VariableNames',{varName})];

        %calculate stage-wise percentages (all of stage, and non-rejected
        %of stage)
        art_pct=calc_percentage(artifact_grids(:,:,strcmp(epochList,currLabel)));
        art_nonrejected_pct=calc_percentage(artifact_grids(:,:,strcmp(epochList,currLabel) & ~reject_grid(1,:)));

        %art_nonrejected_pct=calc_percentage(artifact_grids(:,:,strcmp(scoring_artifact_level.epochs(~reject_grid(1,:)),currLabel)));

        varName=strjoin({'artifact' 'data', currLabel,'perc'},'_');
        varName_2=strjoin({'artifact' 'nonrejected', currLabel,'perc'},'_');

        summaryTable=[summaryTable array2table([art_pct art_nonrejected_pct],'RowNames',artifact_labels,'VariableNames',{varName varName_2})];

    end

    %copy scoring and set excluded to match the reject_mat
    scoring_artifact_level_exclude=scoring_artifact_level;
    scoring_artifact_level_exclude.excluded=reject_grid(1,1:length(scoring_artifact_level_exclude.excluded));

    %scorings
    cfg_artifacts.scorings.scoring_epoch_level=scoring;
    cfg_artifacts.scorings.scoring_artifact_level=scoring_artifact_level;
    cfg_artifacts.scorings.scoring_artifact_level_exclude=scoring_artifact_level_exclude;

end



cfg_artifacts.artifacts.summary=summaryTable;

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
