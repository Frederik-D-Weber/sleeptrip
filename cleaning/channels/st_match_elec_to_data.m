function elec_new=st_match_elec_to_data(cfg)

% ST_MATCH_ELEC_TO_DATA generates an elec configuration matching the
% channels in data
% 
% Use as:
%     elec=st_match_elec_to_data(cfg)
% 
% Required configuration parameters:
%     cfg.elec      = structure, containing ORIGINAL channel information (FieldTrip format)
%     cfg.data      = structure, containing data (FieldTrip format)
%
% Output:
%     elec_new          = structure, containing channel information matching data (FieldTrip format)
%
% See also FT_READ_SENS

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
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

%-----core function start---
%---input checks and defaults----
ft_checkconfig(cfg,'required',{'elec','data'});

fprintf([functionname ' function initialized\n']);

elec_ori=cfg.elec;
data=cfg.data;

%channel labels according to data and (larger) elec cfg
labels_data=data.label;
labels_elec=elec_ori.label;

%all labels in DATA should exist in ELEC_ORI
[existsInElec, selectInds]=ismember(labels_data,labels_elec);
if ~all(existsInElec)
    elec_new=elec_ori;
    fprintf(['WARNING: the following channel labels (from data) could not be found in elec:\n' ...
        '%s\nnot adjusting elec configuration\n'],strjoin(labels_data(~existsInElec)))
    return
else
    elec_new=[];
    num_elec_full=size(labels_elec,1);

    elecFields=fieldnames(elec_ori);
    for field_i=1:length(elecFields)
        fieldName=elecFields{field_i};
        fieldValue=elec_ori.(fieldName);
        if size(fieldValue,1)==num_elec_full
            elec_new.(fieldName)=fieldValue(selectInds,:);
        else
            elec_new.(fieldName)=fieldValue;
        end

    end

end
%---core function end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)