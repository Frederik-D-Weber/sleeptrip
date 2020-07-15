function st = sleepStage2str(st)
%convert sleep stages to abstract level 1

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
    case {'W' 'w' 'Wake' 'wake' 'WAKE' '0'}
        st = 'W';
    case {'S1' 'N1' 'stage 1' 'Stage 1' 'Stage1' 'STAGE 1' 'STAGE1' 'S 1' 'Stadium 1' 'STADIUM 1' 'STADIUM1' '1'}
        st = 'S';
    case {'S2' 'N2' 'stage 2' 'Stage 2' 'Stage2' 'STAGE 2' 'STAGE2' 'S 2' 'Stadium 2' 'STADIUM 2' 'STADIUM2' '2' }
        st = 'S';
    case {'S3' 'N3' 'stage 3' 'Stage 3' 'Stage3' 'STAGE 3' 'STAGE3' 'S 3' 'Stadium 3' 'STADIUM 3' 'STADIUM3' '3' 'SWS'}
        st = 'S';
    case {'S4' 'N4' 'stage 4' 'Stage 4' 'Stage4' 'STAGE 4' 'STAGE4' 'S 4' 'Stadium 4' 'STADIUM 4' 'STADIUM4' '4' 'SWS4'}
        st = 'S';
    case {'REM' 'R' 'r' 'Rem' 'rem' 'Stage 5' 'Stage5' 'STAGE 5' 'STAGE5' 'S 5' 'Stadium 5' 'STADIUM 5' 'STADIUM5' '5'}
        st = 'S';
    case {'MT' 'mt' 'movement' 'Movement' 'Movement Time' 'MovementTime' '8'}
        st = 'W';
    case {'A' 'a' 'Artifact' 'Artefact' 'artifact' 'artefact' 'Artf' 'Artif.'}
        st = 'A';
    otherwise
        st = '?';
end
end