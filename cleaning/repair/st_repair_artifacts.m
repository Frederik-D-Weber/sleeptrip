function data=st_repair_artifacts(cfg_artifacts,data)

% ST_REPAIR_ARTIFACTS interpolates channels segment-wise, as specified by
% the repair grid created by ST_PROCESS_DETECTOR_RESULTS
%
% Use as:
%     data=st_repair_artifacts(cfg,data)
%
% Required configuration parameters:
%     cfg.grid      = structure containing segment-based artifact grids (including repair grid)
%     data          =  data structure to be cleaned
%
% Output:
%     data =  segment-wise interpolated data
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


repair_grid=cfg_grid.repair_grid;
[numChan,numSegments]=size(repair_grid);

elec=cfg_artifacts.elec;

segmentLengthSample=round(cfg_grid.segment_length*data.fsample);

repairSegment_inds=find(any(repair_grid,1));

fprintf('repairing %i of %i segments\n',length(repairSegment_inds),numSegments)

for segment_i=repairSegment_inds
    
    %reimplementation of Fieldtrip's channel interpolation
    badInds=find(repair_grid(:,segment_i)==1);
    goodInds=find(repair_grid(:,segment_i)==0);
    
    repair  = eye(numChan, numChan); %set diagonals to 1
    
    for bad_i=1:length(badInds)
        repair(badInds(bad_i),badInds(bad_i))=0; %set diagonal of bad chan to 0
        distance = sqrt(sum((elec.chanpos(goodInds, :) - repmat(elec.chanpos(badInds(bad_i), :), length(goodInds), 1)).^2, 2));
        
        repair(badInds(bad_i), goodInds) = (1./distance);
        repair(badInds(bad_i), goodInds) = repair(badInds(bad_i), goodInds) ./ sum(repair(badInds(bad_i), goodInds));
    end
    
    %get segment start/ends
    segment_start_sample=(segment_i-1)*segmentLengthSample+1;
    segment_end_sample=segment_i*segmentLengthSample;
    
    %check we're inside the data range
    if segment_end_sample>data.sampleinfo(2)
        segment_end_sample=data.sampleinfo(2);
    end
    
    %replace original data with repair*original (matrix multiplication)
    data.trial{1}(:,segment_start_sample:segment_end_sample)=repair*data.trial{1}(:,segment_start_sample:segment_end_sample);
end

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)