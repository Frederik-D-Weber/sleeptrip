function st = sleepStage2str_alt(st)
%convert sleep stages to abstract level 2

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
switch st
    case {'S3' 'N3' 'stage 3' 'Stage 3' 'Stage3' 'STAGE 3' 'STAGE3' 'S 3' 'Stadium 3' 'STADIUM 3' 'STADIUM3' '3' 'SWS'}
        st = 'SWS';
    case {'S4' 'N4' 'stage 4' 'Stage 4' 'Stage4' 'STAGE 4' 'STAGE4' 'S 4' 'Stadium 4' 'STADIUM 4' 'STADIUM4' '4' 'SWS4'}
        st = 'SWS';
    otherwise
        st = sleepStage2str(st);
end
end