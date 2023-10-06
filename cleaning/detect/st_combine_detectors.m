function detector_set=st_combine_detectors(all_cfg)

% ST_COMBINE_DETECTORS combines indvidual detector configurations into
% a "detector set" (=collection of artifact detectors),
% which can be supplied to st_run_detector_set.
%
% Use as:
%     detector_set=st_combine_detectors(all_cfg)
%
% Required input:
%     all_cfg      = cell array of individual artifact detector cfgs
%
% Output:
%     detector_set  = structure/cfg, containing details of individual
%     detectors
%
% See also ST_GET_DEFAULT_DETECTOR_SET

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


%------start core function
%---input checks and defaults----
inputOk=isa(all_cfg, 'cell') && all(cellfun(@isstruct, all_cfg(:)));
if ~inputOk
    ft_error('Input should be a cell array of structures.\n');
end

fprintf([functionname ' function initialized\n']);

%put all individual detectors in an output cfg
detector_set=[];
detector_set.number=length(all_cfg);
detector_set.label=cellfun(@(X) X.label,all_cfg,'UniformOutput',false);
detector_set.detectors=all_cfg;
%---end core function


% do the general cleanup and bookkeeping at the end of the function

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)