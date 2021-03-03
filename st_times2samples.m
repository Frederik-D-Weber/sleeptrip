function samples = st_times2samples(data, times_seconds, rmnanrows)

% ST_TIMES2SAMPLES converts time values from seconds to the samples 
% that match the given data structure time. The data needs to be a
% continous data
%
% Use as
%   samples = st_convertsamples(data, times_seconds)
%   samples = st_convertsamples(data, times_seconds, 'rmnanrows')
%
%
%
% See also FT_DEFINETRIAL, FT_REDEFINETRIAL, FT_PREPROCESSING

% Copyright (C) 2019-, Frederik D. Weber
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

% data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'hassampleinfo', 'yes');
% if isfield(data, 'trial') && isfield(data, 'time')
%     % check if data structure is likely continous data
%     if ((numel(data.trial) ~= 1) || ~all(size(data.sampleinfo) == [1 2])) && ~(nargout > 1)
%         ft_error('data structure does not look like continous data and has more than one trial')
%     end
% end

if iscell(data.time)
    tsref = data.time{1};
else
    tsref = data.time;    
end

samples = interp1(tsref, 1:numel(tsref), times_seconds, 'nearest');

no_match_count = sum(sum(isnan(samples)));
if no_match_count>0
    ft_warning('%d times_seconds could not be matched.\n please check if they fall within the data boundary.',no_match_count)
    if (nargin>2) && strcmp(rmnanrows,'rmnanrows')
        samples = samples(sum(isnan(samples),2)==0,:);
    end
end

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end