function scoring_artifact_level=st_scoring_artifact_level(cfg)

% ST_SCORING_ARTIFACT_LEVEL converts an original scoring (typically with 30 s epochs)
% to a scoring with epoch lengths matching the artifact segment length (typically 5 s),
% and marking epochs for exclusion based on artifact rejection.
% This artifact-level scoring can then be used to select clean data of particular stages.
%
% Use as:
%     scoring_artifact_level=st_scoring_artifact_level(cfg)
%
% Required configuration parameters:
%       cfg.scoring  = a scoring structure as defined in ST_READ_SCORING
%       cfg.cfg_artifacts  = a cfg_artifacts structure as returned by ST_PROCESS_DETECTOR_RESULTS
%
% Optional configuration parameters:
%       cfg.userejection = 'yes' or 'no'. Whether or not to ignore include rejection information (default: 'yes').
%
%Output:
%     scoring_artifact_level = a scoring structure as defined in ST_READ_SCORING
%
% See also ST_READ_SCORING ST_PROCESS_DETECTOR_RESULTS ST_SELECT_DATA

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
ft_checkconfig(cfg,'required',{'scoring','cfg_artifacts'});
cfg.userejection=ft_getopt(cfg,'userejection','yes');


scoring = cfg.scoring;
cfg_artifacts=cfg.cfg_artifacts;

%---level deeper:
% input checks and defaults----
ft_checkconfig(cfg_artifacts.artifacts,'required',{'grid'});
cfg_grid=cfg_artifacts.artifacts.grid;



%get basics
segment_length=cfg_grid.segment_length;
numArtifactSegment=cfg_grid.segment_number;

%rescore to match artifact segments
tmp_cfg=[];
tmp_cfg.new_epochlength_seconds=segment_length;
scoring_artifact_level = st_scoring_change_epochlength(tmp_cfg, scoring);

numScoringSegment=length(scoring_artifact_level.epochs);

fprintf('Changing epoch length to match artifact segment length: %i epochs versus %i artifact segments\n',numScoringSegment,numArtifactSegment)

segmentDisagree=abs(numScoringSegment-numArtifactSegment);
minutesDisagree=segmentDisagree*segment_length/60;

if minutesDisagree>=1
    ft_warning('total duration of epochs and artifact segments differ by %.1f minutes: consider checking',minutesDisagree)
end

%artifact-level scoring may be longer or shorter than artifact segments:
if numArtifactSegment<numScoringSegment % trim scoring if needed

    fprintf('Removing %i final epochs from scoring\n',numScoringSegment-numArtifactSegment)
    scoring_artifact_level.epochs=scoring_artifact_level.epochs(1:numArtifactSegment);
    scoring_artifact_level.excluded=scoring_artifact_level.excluded(1:numArtifactSegment);
    scoring_artifact_level.prob=scoring_artifact_level.prob(:,1:numArtifactSegment);
    scoring_artifact_level.numbers=scoring_artifact_level.numbers(1:numArtifactSegment);
elseif numArtifactSegment>numScoringSegment %extend scoring by repeating last epoch info

    fprintf('Adding %i epochs to end of scoring by repeating final epoch\n',numArtifactSegment-numScoringSegment)
    scoring_artifact_level.epochs(end+1:numArtifactSegment)=scoring_artifact_level.epochs(end);
    scoring_artifact_level.excluded(end+1:numArtifactSegment)= scoring_artifact_level.excluded(end);
    scoring_artifact_level.prob(:,end+1:numArtifactSegment)=repmat(scoring_artifact_level.prob(:,end),[1 numArtifactSegment-numScoringSegment]);
    scoring_artifact_level.numbers(end+1:numArtifactSegment)= scoring_artifact_level.numbers(end);
end

%mark epochs for exclusion
if istrue(cfg.userejection)
    scoring_artifact_level.excluded=cfg_grid.reject_grid(1,:);

end

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)